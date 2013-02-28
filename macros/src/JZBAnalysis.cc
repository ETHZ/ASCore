#include "helper/Utilities.hh"
#include "JZBAnalysis.hh"
#include "TH1.h"
#include <time.h>
#include <TRandom.h>
#include "TF1.h"
#include "TTree.h"
#include <cstdlib>
using namespace std;


#define jMax 30  // do not touch this
#define rMax 30
#define Zmax 30

enum METTYPE { mettype_min, RAW = mettype_min, T1PFMET, TCMET, MUJESCORRMET, PFMET, SUMET, PFRECOILMET, RECOILMET, mettype_max };
enum JZBTYPE { jzbtype_min, TYPEONECORRPFMETJZB = jzbtype_min, PFJZB, RECOILJZB, PFRECOILJZB, TCJZB, jzbtype_max };

string sjzbversion="$Revision: 1.70.2.119 $";
string sjzbinfo="";
TRandom3 *r;
TF1 *L5corr_bJ;
TF1 *L5corr_qJ;
TF1 *L5corr_cJ;
TF1 *L5corr_gJ;


float firstLeptonPtCut  = 10.0;
float secondLeptonPtCut = 10.0;

/*
$Id: JZBAnalysis.cc,v 1.70.2.119 2013/02/21 10:37:22 buchmann Exp $
*/


Double_t GausRandom(Double_t mu, Double_t sigma) { 
  return gRandom->Gaus(mu,sigma);   //real deal
  //return mu;//debugging : no smearing
}

class nanoEvent
{
public:
  nanoEvent();
  void reset();

  float mll; // di-lepton system
  float minc;
  float pt;
  float phi;
  float eta;
  float E;
  bool is_data;

  float l1l2dR;
  float pt1; // leading leptons
  float pt2;
  float iso1;
  float iso2;
  bool  isConv1; // Photon conversion flag
  bool  isConv2;
  bool softMuon;
  bool softMuonMC;

  int NgenLeps;
  int NgenZs;
  
  float genPt1; // leading legenPtons
  float genPt2;
  float genDRll;
  float genPhi1;
  float genPhi2;
  int   genId1;
  int   genId2;
  int   genMID1gen;
  int   genMID2gen;
  int   genMID1;
  int   genMID2;
  int   genGMID1gen;
  int   genGMID2gen;
  int   genGMID1;
  int   genGMID2;
  float genEta1; // leading legenPtons
  float genEta2;
  float genMET;
  float genZPt;    // True Z Pt
  float genMll;
  float genRecoil;
  float genJZB;
  int   genNjets;
  int   genNleptons;
  float genRecoilSel;
  float genPt1Sel; // Selected leptons
  float genPt2Sel;
  float genEta1Sel;
  float genEta2Sel;
  int   genId1Sel;
  int   genId2Sel;
  float genZPtSel; // Z candidate from selected leptons
  float genMllSel;
  float genJZBSel;
  float eta1; // leading leptons
  float eta2;
  float phi1;
  float phi2;
  float dphi;
  float dphiZpfMet;
  float dphigenZgenMet;
  float dphiZs1;
  float dphiZs2;
  float dphiMet1;
  float dphiMet2;
  float dphitcMet1;
  float dphitcMet2;
  float dphipft1Met1;
  float dphipft1Met2;
  float dphipfRecoilMet1;
  float dphipfRecoilMet2;
  
  float pfMET;
  float t1pfMET;

  bool ElCInfoIsGsfCtfCons;
  bool ElCInfoIsGsfCtfScPixCons;
  bool ElCInfoIsGsfScPixCons;

  int id1;
  int id2;
  int ch1;
  int ch2;
  int chid1; // old id (kostas convention)
  int chid2;

  int process;

  int leptonNum; // store all leptons (reduntant for the 2 leptons that form the Z)
  float leptonPt[jMax]; 
  float leptonEta[jMax];
  float leptonPhi[jMax];
  int leptonId[jMax];
  int leptonCharge[jMax];

  int leptonPairNum;
  int leptonPairId[jMax];
  float leptonPairMass[jMax];
  float leptonPairDphi[jMax];


  int pfJetNum;
  float pfJetPt[jMax];
  float pfJetEta[jMax];
  float pfJetPhi[jMax];
  bool  pfJetID[jMax];
  float pfJetScale[jMax];
  float pfJetScaleUnc[jMax];
  float pfJetDphiMet[jMax];
  float pfJetDphiZ[jMax];
  float pfBJetDphiZ[jMax];
  float CorrectionRatio[jMax];
  float pfHT;
  float pfGoodHT;
  float pfTightHT;
  
  float metUncertainty;
  float type1metUncertainty;
  
  int pfJetGoodNum30;
  int pfJetGoodNum30Fwd;
  int pfJetGoodNumID;
  int pfJetGoodNump1sigma;
  int pfJetGoodNumn1sigma;
  int pfJetGoodNum40p1sigma;
  int pfJetGoodNum40n1sigma;
  int pfJetGoodNum50p1sigma;
  int pfJetGoodNum50n1sigma;
  int pfJetGoodNum40CHS;
  float pfJetGoodPt[jMax];
  float pfJetGoodEta[jMax];
  float pfJetGoodPhi[jMax];
  float pfJetGoodE[jMax];
  float pfJetGoodMl[jMax];
  float pfJetGoodM[jMax];
  float pfJetGoodPtl[jMax];
  float pfJetGoodID[jMax];
  int   pfJetGoodTracks[jMax];
  int   pfJetGoodTracksN[jMax];
  float bTagProbCSVBP[jMax];
  float bTagProbCSVMVA[jMax];
  
  int pfJetGoodNumBtag30;
  int pfJetGoodNumBtag40;
  int pfJetGoodNumBtag30_Tight;
  int pfJetGoodNumBtag40_Tight;
  int pfJetGoodNumBtag30_Loose;
  int pfJetGoodNumBtag40_Loose;
  int pfJetGoodNumIDBtag;
  float pfJetGoodPtBtag[jMax];
  float pfJetGoodEtaBtag[jMax];
  float pfJetGoodPhiBtag[jMax];
  float pfJetGoodEBtag[jMax];
  float pfJetGoodMlBtag[jMax];
  float pfJetGoodMBtag[jMax];
  float pfJetGoodPtlBtag[jMax];
  bool  pfJetGoodIDBtag[jMax];
  int   pfJetGoodTracksBtag[jMax];
  int   pfJetGoodTracksNBtag[jMax];

  int pfJetGoodNum40;
  int pfJetGoodNum40Fwd;
  int pfJetGoodNum50;
  int pfJetGoodNum50Fwd;
  int pfJetGoodNum60;
  int pfJetGoodNum60Fwd;

  int EventFlavor;
  bool EventZToTaus;
  
  float met[mettype_max];
  float metPhi[mettype_max];
  float dphiMetLep[mettype_max];
  float dphiMetJet[mettype_max];
  float dphiMetSumJetPt[mettype_max];
  float metPerp[mettype_max];
  float metPar[mettype_max];
  ULong64_t eventNum;
  int runNum;
  int lumi;
  int goodVtx;
  int numVtx;
  float totEvents; // tot events processed by the ntuple producer (job submission efficiency), no need to keep this as int, better as float
  int badJet;

  float jzb[jzbtype_max];
  float sjzb[jzbtype_max]; // smeared JZB
  float dphi_sumJetVSZ[jzbtype_max];
  float sumJetPt[jzbtype_max];

  float weight;
  float weightEffDown;
  float weightEffUp;
  float Efficiencyweightonly;

  
  int NPdfs;
  float pdfW[100];
  float pdfWsum;
  float PUweight,PUweightUP,PUweightDOWN;
  bool passed_triggers;
  int trigger_bit;
  bool passed_filters;
  int filter_bit;
  float mGlu;
  float mChi;
  float mLSP;
  float xSMS;
  float xbarSMS;
  float mGMSBGlu;
  float mGMSBChi;
  float mGMSBLSP;
  
  
  //gen information
  int nZ; // number of generator Z's in the process
  int SourceOfZ[Zmax];//mother particle of the (first Zmax) Z's
  int DecayCode; //decay code: 100*h + l, where h = number of hadronically decaying Z's, l = number of leptonically decaying Z's (e.g. 102 = 1 had. Z, 2 lep. Z's)
  float realx; // this is the "x" we measure (for scans)
  float imposedx; // this is the "x" we imposed.
  float pureGeneratorZpt;
  float pureGeneratorZM;
  float pureGeneratorZphi;
  float pureGeneratorZeta;
  float pureGeneratorJZB;
  float pureGeneratorMet;
  float pureGeneratorMetPhi;
  float pureGeneratorSumJetPt;
  float pureGeneratorSumJetEta;
  float pureGeneratorSumJetPhi;
  float pure2ndGeneratorJZB;
  float pure2ndGeneratorZpt;
  int nLSPs;
  float angleLSPLSP;
  float angleLSPLSP2d;
  float angleLSPZ;
  float angleLSPZ2d;
  float angleChi2Z2d;
  float angleChi2Z;
  float dphiSumLSPgenMET;
  float SumLSPEta;
  float SumLSPPhi;
  float absvalSumLSP;
  int LSPPromptnessLevel[2];
  int ZPromptnessLevel[2];
  float LSP1pt;
  float LSP2pt;
  int LSP1Mo;
  int LSP2Mo;
  float LSP1Mopt;
  float LSP2Mopt;
  
  //Z+b variables
  
  bool HasSoftLepton;
  float SoftLeptonPt;
  
  float mpf;
  float fake_mpf;
  
  float ZbCHS3010_alpha;
  float ZbCHS3010_alphaUp;
  float ZbCHS3010_alphaDown;
  float ZbCHS3010_alphaL5;
  float ZbCHS3010_L5corr;
  
  float ZbCHS1010_alpha;
  float ZbCHS1010_alphaUp;
  float ZbCHS1010_alphaDown;
  
  float ZbCHS3010_bTagProbCSVBP[jMax]; 
  int ZbCHS3010_pfJetGoodNumBtag; 
  float ZbCHS3010_pfJetGoodEta[jMax]; 
  int ZbCHS3010_pfJetGoodNum;
  float ZbCHS3010_pfJetDphiZ[jMax];
  float ZbCHS3010_pfJetGoodPt[jMax];
  float ZbCHS3010_pfJetSum;
  float ZbCHS3010_pfBJetDphiZ[jMax];

  float ZbCHS1010_bTagProbCSVBP[jMax]; 
  int ZbCHS1010_pfJetGoodNumBtag; 
  float ZbCHS1010_pfJetGoodEta[jMax]; 
  int ZbCHS1010_pfJetGoodNum;
  float ZbCHS1010_pfJetDphiZ[jMax];
  float ZbCHS1010_pfJetGoodPt[jMax];
  float ZbCHS1010_pfJetSum;
  float ZbCHS1010_pfBJetDphiZ[jMax];

  float ZbCHS3010_BTagWgtT;
  float ZbCHS3010_BTagWgtTUp;
  float ZbCHS3010_BTagWgtTDown;

  float ZbCHS1010_BTagWgtT;
  float ZbCHS1010_BTagWgtTUp;
  float ZbCHS1010_BTagWgtTDown;
  
  float ZbCHS3010_BTagWgtM;
  float ZbCHS3010_BTagWgtMUp;
  float ZbCHS3010_BTagWgtMDown;
  
  float ZbCHS1010_BTagWgtM;
  float ZbCHS1010_BTagWgtMUp;
  float ZbCHS1010_BTagWgtMDown;
  
  float ZbCHS3010_BTagWgtL;
  float ZbCHS3010_BTagWgtLUp;
  float ZbCHS3010_BTagWgtLDown;
  
  float ZbCHS1010_BTagWgtL;
  float ZbCHS1010_BTagWgtLUp;
  float ZbCHS1010_BTagWgtLDown;
  
  bool ZbCHS3010_LeadingJetIsPu;
  bool ZbCHS1010_LeadingJetIsPu;
  
  bool ZbCHS3010_SubLeadingJetIsPu;
  bool ZbCHS1010_SubLeadingJetIsPu;
  
  
  float fact;
  
 //Trilepton variables
  float tri_pt1;
  float tri_pt2;
  float tri_pt3;
  float tri_eta1;
  float tri_eta2;
  float tri_eta3;
  bool tri_MatchFound;
  int tri_id1;
  int tri_id2;
  int tri_id3;
  int tri_ch1;
  int tri_ch2;
  int tri_ch3;
  float tri_mlll;
  float tri_badmll;
  float tri_mll;
  float tri_submll;
  float tri_badsubmll;
  float tri_mT;
  float tri_badmT;
  int tri_genMID1;
  int tri_genMID2;
  int tri_genMID3;
  int tri_badgenMID1;
  int tri_badgenMID2;
  int tri_badgenMID3;

  float tri_dR12;
  float tri_dR13;
  float tri_dR23;
  int tri_index1;
  int tri_index2;
  int tri_index3;
  int tri_badindex1;
  int tri_badindex2;
  int tri_badindex3;
  int tri_badid1;
  int tri_badid2;
  int tri_badid3;
  bool tri_GoodZMatch;
  bool tri_GoodWMatch;
  
  float gentri_pt1;
  float gentri_pt2;
  float gentri_pt3;
  float gentri_eta1;
  float gentri_eta2;
  float gentri_eta3;
  float gentri_id1;
  float gentri_id2;
  float gentri_id3;
  float gentri_ch1;
  float gentri_ch2;
  float gentri_ch3;
  float gentri_mlll;
  float gentri_mll;
  float gentri_badmll;
  float gentri_submll;
  float gentri_badsubmll;
  float gentri_mT;
  float gentri_badmT;
  int gentri_genMID1;
  int gentri_genMID2;
  int gentri_genMID3;
  int gentri_badgenMID1;
  int gentri_badgenMID2;
  int gentri_badgenMID3;
  float gentri_dR12;
  float gentri_dR13;
  float gentri_dR23;
  bool gentri_GoodZMatch;
  bool gentri_GoodWMatch;
  
  float pgentri_pt1;
  float pgentri_pt2;
  float pgentri_pt3;
  float pgentri_eta1;
  float pgentri_eta2;
  float pgentri_eta3;
  float pgentri_id1;
  float pgentri_id2;
  float pgentri_id3;
  float pgentri_ch1;
  float pgentri_ch2;
  float pgentri_ch3;
  float pgentri_mlll;
  float pgentri_mll;
  float pgentri_submll;
  float pgentri_mT;
  float pgentri_dR12;
  float pgentri_dR13;
  float pgentri_dR23;
  
  float dRmerged;
  int merged;
  float mergedPt;
  float mergedEta;
  float mergedPhoHCalIso2012ConeDR03;
  float mergedPhoSigmaIetaIeta;
  float mergedpfchargedhadiso;
  float mergedpfneutralhadiso;
  float mergedpfphotoniso;

};

nanoEvent::nanoEvent(){};
void nanoEvent::reset()
{

  mll=0; // di-lepton system
  minc=0;
  pt=0;
  phi=0;
  eta=0;
  E=0;

  is_data=false;
  NPdfs=0;
  pdfWsum=0;

  process=0;

  l1l2dR=-99.9;
  pt1=0;
  pt2=0;
  iso1=0;
  iso2=0;
  isConv1 = false;
  isConv2 = false;
  softMuon = false;
  softMuonMC = false;
  
  NgenZs=0;
  NgenLeps=0;

  genPt1=0;
  genPt2=0;
  genDRll=0;
  genPhi1=0;
  genPhi2=0;
  genEta1=0;
  genEta2=0;
  genId1=0;
  genId2=0;
  genMID1gen=0;
  genMID2gen=0;
  genMID1=0;
  genMID2=0;
  genGMID1gen=0;
  genGMID2gen=0;
  genGMID1=0;
  genGMID2=0;
  genMET=0;
  genZPt=0;
  genMll=0;
  genRecoil=0;
  genJZB = 0;
  genNjets = 0;
  genNleptons = 0;
  genPt1Sel=0;
  genPt2Sel=0;
  genEta1Sel=0;
  genEta2Sel=0;
  genId1Sel=0;
  genId2Sel=0;
  genZPtSel=0;
  genMllSel=0;
  genRecoilSel=0;
  genJZBSel = 0;
  passed_triggers=0;
  trigger_bit = 0;
  passed_filters=0;
  filter_bit = 0;
  
  eta1=0; // leading leptons
  eta2=0;
  phi1=0;
  phi2=0;
  dphiZpfMet=0;
  dphigenZgenMet=0;
  dphiZs1=0;
  dphiZs2=0;
  dphiMet1=0;
  dphiMet2=0;
  dphitcMet1=0;
  dphitcMet2=0;
  dphipft1Met1=0;
  dphipft1Met1=0;
  dphipfRecoilMet1=0;
  dphipfRecoilMet2=0;
  dphi=0;
  ElCInfoIsGsfCtfCons=false;
  ElCInfoIsGsfCtfScPixCons=false;
  ElCInfoIsGsfScPixCons=false;
  id1=-9;
  id2=-9;
  ch1=-9;
  ch2=-9;
  chid1=0;
  chid2=0;

  for(int i=0;i<100;i++) pdfW[i]=1.0;

  for(int jCounter=0;jCounter<jMax;jCounter++){
    leptonPt[jCounter]=0; 
    leptonEta[jCounter]=0;
    leptonPhi[jCounter]=0;
    leptonId[jCounter]=0;
    leptonCharge[jCounter]=0;
  }
  leptonNum=0;

  for(int jCounter=0;jCounter<jMax;jCounter++){
    leptonPairMass[jCounter]=0;
    leptonPairDphi[jCounter]=0;
    leptonPairId[jCounter]=0;
  } 
  leptonPairNum=0;
 
  for(int metCounter=int(mettype_min);metCounter<int(mettype_max);metCounter++){
    met[metCounter]=0;
    metPhi[metCounter]=0;
    dphiMetLep[metCounter]=0;
    dphiMetJet[metCounter]=0;
    dphiMetSumJetPt[metCounter]=0;
    metPerp[metCounter]=0;
    metPar[metCounter]=0;
   
  }

  for(int jCounter=0;jCounter<jMax;jCounter++){
    pfJetPt[jCounter]=0;
    pfJetEta[jCounter]=0;
    pfJetPhi[jCounter]=0;
    pfJetID[jCounter]=0;
    pfJetScale[jCounter]=0;
    pfJetScaleUnc[jCounter]=0;
    pfJetDphiMet[jCounter]=0;
    pfJetDphiZ[jCounter]=0;
    pfBJetDphiZ[jCounter]=0;
    CorrectionRatio[jMax]=0;
  }
  pfJetNum=0;
  pfHT=0;
  pfGoodHT=0;
  pfTightHT=0;

  for(int jCounter=0;jCounter<jMax;jCounter++){
    pfJetGoodPt[jCounter]=0;
    pfJetGoodEta[jCounter]=0;
    pfJetGoodPhi[jCounter]=0;
    pfJetGoodE[jCounter]=0;
    pfJetGoodM[jCounter]=0;
    pfJetGoodMl[jCounter]=0;
    pfJetGoodPtl[jCounter]=0;
    pfJetGoodID[jCounter]=0;
    pfJetGoodTracks[jCounter]=0;
    pfJetGoodTracksBtag[jCounter]=0;
    pfJetGoodTracksN[jCounter]=0;
    pfJetGoodTracksNBtag[jCounter]=0;
    pfJetGoodEtaBtag[jCounter]=0;
    pfJetGoodPhiBtag[jCounter]=0;
    pfJetGoodEBtag[jCounter]=0;
    pfJetGoodMlBtag[jCounter]=0;
    pfJetGoodMBtag[jCounter]=0;
    pfJetGoodPtlBtag[jCounter]=0;
    pfJetGoodIDBtag[jCounter]=0;
    bTagProbCSVBP[jCounter] = 0;
    bTagProbCSVMVA[jCounter] = 0;
  }

  EventFlavor=0;
  EventZToTaus=false;
  
  metUncertainty=0.0;
  type1metUncertainty=0.0;

  
  pfJetGoodNum30=0;
  pfJetGoodNumID=0;
  pfJetGoodNumBtag30=0;
  pfJetGoodNumBtag40=0;
  pfJetGoodNumBtag30_Tight=0;
  pfJetGoodNumBtag40_Tight=0;
  pfJetGoodNumBtag30_Loose=0;
  pfJetGoodNumBtag40_Loose=0;
  pfJetGoodNumIDBtag=0;
  pfJetGoodNump1sigma=0;
  pfJetGoodNumn1sigma=0;
  pfJetGoodNum40p1sigma=0;
  pfJetGoodNum40n1sigma=0;
  pfJetGoodNum50p1sigma=0;
  pfJetGoodNum50n1sigma=0;

  pfJetGoodNum30Fwd=0;
  
  pfJetGoodNum40CHS=0;
  pfJetGoodNum40=0;
  pfJetGoodNum40Fwd=0;
  pfJetGoodNum50=0;
  pfJetGoodNum50Fwd=0;
  pfJetGoodNum60=0;
  pfJetGoodNum60Fwd=0;

  eventNum=0;
  runNum=0;
  lumi=0;
  goodVtx=0;
  numVtx=0;
  badJet=0;
  totEvents=0;

  for(int rCounter=int(jzbtype_min);rCounter<int(jzbtype_max);rCounter++){
    jzb[rCounter]=0;
    sjzb[rCounter]=0;
    dphi_sumJetVSZ[rCounter]=0;
    sumJetPt[rCounter]=0;
  }

  weight = 1.0;
  PUweight = 1.0;
  PUweightUP = 1.0;
  PUweightDOWN = 1.0;
  Efficiencyweightonly = 1.0;
  weightEffDown = 1.0;
  weightEffUp = 1.0;

  ZbCHS3010_BTagWgtT = 1.0 ;
  ZbCHS1010_BTagWgtT = 1.0 ;
  
  ZbCHS3010_BTagWgtTUp = 1.0 ;
  ZbCHS1010_BTagWgtTUp = 1.0 ;
  
  ZbCHS3010_BTagWgtTDown = 1.0 ;
  ZbCHS1010_BTagWgtTDown = 1.0 ;

  ZbCHS3010_BTagWgtM = 1.0 ;
  ZbCHS1010_BTagWgtM = 1.0 ;
  
  ZbCHS3010_BTagWgtMUp = 1.0 ;
  ZbCHS1010_BTagWgtMUp = 1.0 ;
  
  ZbCHS3010_BTagWgtMDown = 1.0 ;
  ZbCHS1010_BTagWgtMDown = 1.0 ;
  
  ZbCHS3010_BTagWgtL = 1.0 ;
  ZbCHS1010_BTagWgtL = 1.0 ;
  
  ZbCHS3010_BTagWgtLUp = 1.0 ;
  ZbCHS1010_BTagWgtLUp = 1.0 ;
  
  ZbCHS3010_BTagWgtLDown = 1.0 ;
  ZbCHS1010_BTagWgtLDown = 1.0 ;
  
  ZbCHS3010_LeadingJetIsPu = false;
  ZbCHS1010_LeadingJetIsPu = false;

  ZbCHS3010_SubLeadingJetIsPu = false;
  ZbCHS1010_SubLeadingJetIsPu = false;

  
  mGlu=0;
  mChi=0;
  mLSP=0;
  xSMS=0;
  xbarSMS=0;
  mGMSBGlu=0;
  mGMSBChi=0;
  mGMSBLSP=0;
  
  // gen info
  nZ=0;
  for(int i=0;i<Zmax;i++) SourceOfZ[i]=0;
  DecayCode=0;
  realx=0;
  pureGeneratorJZB=0;
  pureGeneratorMet=0;
  pureGeneratorMetPhi=0;
  pureGeneratorSumJetPt=0;
  pureGeneratorSumJetEta=0;
  pureGeneratorSumJetPhi=0;

  pure2ndGeneratorJZB=0;
  pure2ndGeneratorZpt=0;
  pureGeneratorZpt=0;
  pureGeneratorZM=0;
  pureGeneratorZeta=0;
  pureGeneratorZphi=0;
  nLSPs=0;
  angleLSPLSP=0;
  angleLSPLSP2d=-5;
  angleLSPZ=0;
  dphiSumLSPgenMET=0;
  SumLSPEta=0;
  SumLSPPhi=0;
  absvalSumLSP=0;
  angleLSPZ2d=-5;
  angleChi2Z2d=-5;
  angleChi2Z=-5;

  LSPPromptnessLevel[0]=-1;
  LSPPromptnessLevel[1]=-1;
  ZPromptnessLevel[0]=-1;
  ZPromptnessLevel[1]=-1;
  LSP1pt=0;
  LSP2pt=0;
  LSP1Mo=0;
  LSP2Mo=0;
  LSP1Mopt=0;
  LSP2Mopt=0;
  
  //Z+b variables
  ZbCHS3010_alpha=0;
  ZbCHS3010_alphaL5=0;
  ZbCHS3010_L5corr=0;
  ZbCHS1010_alpha=0;
  ZbCHS3010_alphaUp=0;
  ZbCHS1010_alphaUp=0;
  ZbCHS3010_alphaDown=0;
  ZbCHS1010_alphaDown=0;
  mpf=0;
  fake_mpf=0;
  
  HasSoftLepton=false;
  SoftLeptonPt=0.;
  
  
  ZbCHS3010_pfJetGoodNumBtag=0;
  ZbCHS3010_pfJetGoodNum=0;
  ZbCHS3010_pfJetSum=0;

  ZbCHS1010_pfJetGoodNumBtag=0;
  ZbCHS1010_pfJetGoodNum=0;
  ZbCHS1010_pfJetSum=0;
  


  for(int i=0;i<jMax;i++) {
    ZbCHS3010_bTagProbCSVBP[i] = 0;
    ZbCHS3010_pfJetGoodEta[i] = 0;
    ZbCHS3010_pfJetDphiZ[i] = 0;
    ZbCHS3010_pfJetGoodPt[i] = 0;
    ZbCHS3010_pfBJetDphiZ[i] = 0;

    ZbCHS1010_bTagProbCSVBP[i] = 0;
    ZbCHS1010_pfJetGoodEta[i] = 0;
    ZbCHS1010_pfJetDphiZ[i] = 0;
    ZbCHS1010_pfJetGoodPt[i] = 0;
    ZbCHS1010_pfBJetDphiZ[i] = 0;
  }
  
  tri_pt1=0;
  fact=1.0;
  tri_pt2=0;
  tri_pt3=0;
  
  tri_eta1=0;
  tri_eta2=0;
  tri_eta3=0;
  
  tri_MatchFound = false;
    
  tri_id1=0;
  tri_id2=0;
  tri_id3=0;
    
  tri_ch1=0;
  tri_ch2=0;
  tri_ch3=0;
    
  tri_mlll=0;
  tri_mll=0;
  tri_submll=0;
  
  tri_badmll=0;
  tri_badsubmll=0;
  tri_badmT=0;
  tri_badgenMID1=0;
  tri_badgenMID2=0;
  tri_badgenMID3=0;
  tri_badindex1=0;
  tri_badindex2=0;
  tri_badindex3=0;
  tri_badid1=0;
  tri_badid2=0;
  tri_badid3=0;
  
    
  tri_mT=0;
  tri_index1=0;
  tri_index1=0;
  tri_index1=0;
  
  tri_genMID1=0;
  tri_genMID2=0;
  tri_genMID3=0;
  tri_dR12=0;
  tri_dR13=0;
  tri_dR23=0;
  
  gentri_pt1=0;
  gentri_pt2=0;
  gentri_pt3=0;
  
  gentri_id1=0;
  gentri_id2=0;
  gentri_id3=0;
    
  gentri_ch1=0;
  gentri_ch2=0;
  gentri_ch3=0;
    
  gentri_mlll=0;
  gentri_mll=0;
  gentri_submll=0;
    
  gentri_mT=0;
  
  gentri_genMID1=0;
  gentri_genMID2=0;
  gentri_genMID3=0;
  gentri_dR12=0;
  gentri_dR13=0;
  gentri_dR23=0;
  
  gentri_badmll=0;
  gentri_badsubmll=0;
  gentri_badmT=0;
  gentri_badgenMID1=0;
  gentri_badgenMID2=0;
  gentri_badgenMID3=0;
  
  gentri_GoodZMatch=false;
  gentri_GoodWMatch=false;
  
  pgentri_pt1=0;
  pgentri_pt2=0;
  pgentri_pt3=0;
  pgentri_eta1=0;
  pgentri_eta2=0;
  pgentri_eta3=0;
  pgentri_id1=0;
  pgentri_id2=0;
  pgentri_id3=0;
  pgentri_ch1=0;
  pgentri_ch2=0;
  pgentri_ch3=0;
  pgentri_mlll=0;
  pgentri_mll=0;
  pgentri_submll=0;
  pgentri_mT=0;
  pgentri_dR12=0;
  pgentri_dR13=0;
  pgentri_dR23=0;

  dRmerged=0.;
  merged=0;
  mergedPt=0;
  mergedEta=0;
  mergedPhoHCalIso2012ConeDR03=0;
  mergedPhoSigmaIetaIeta=0;
  mergedpfchargedhadiso=0;
  mergedpfneutralhadiso=0;
  mergedpfphotoniso=0;
  
}


TTree *myTree;
TTree *myInfo;
TH1F *weight_histo;

nanoEvent nEvent;

bool IsLepton(const int pdgid) {
  int tid = abs(pdgid);
  if(tid==11||tid==13||tid==15) return true;
  else return false;
}

bool IsBMeson(const int pdgid) {
  int MesonIDs[99]={511,521,10511,10521,513,523,10513,10523,20513,20523,515,525,531,10531,533,10533,20533,535,541,10541,543,10543,20543,545,5101,5103,5201,5203,5301,5303,5401,5403,5503,10113,10213,551,10551,100551,110551,200551,210551,553,10553,20553,30553,100553,110553,120553,130553,200553,210553,220553,300553,9000553,9010553,555,10555,20555,100555,110555,120555,200555,557,100557,5122,5112,5212,5222,5114,5214,5224,5132,5232,5312,5322,5314,5324,5332,5334,5142,5242,5412,5422,5414,5424,5342,5432,5434,5442,5444,5512,5522,5514,5524,5532,5534,5542,5544,5554};
  
  for(int i=0;i<99;i++) {
    if(pdgid==MesonIDs[i]) return true;
  }
  
  return false;
}

void JZBAnalysis::IsParticleFromB(const int index) {
  if(index<0) return;
  if(fTR->genInfoId[index]==23) return;
  if(IsBMeson(fTR->genInfoId[index])) {
    nEvent.HasSoftLepton=true;
    return;
  } else {
    int motherindex=fTR->genInfoMo1[index];
    int motherindex2=fTR->genInfoMo2[index];
    if(motherindex>=0) IsParticleFromB(motherindex);
    if(motherindex2>=0) IsParticleFromB(motherindex2);
  }
}

void JZBAnalysis::ContainsSoftLepton() {
  if(fTR->nGenParticles<1) return;
  for(int i=0;i<fTR->nGenParticles;i++) {
    int thisParticleId = fTR->genInfoId[i];
    if(IsLepton(thisParticleId)) {
      int motherindex=fTR->genInfoMo1[i];
      int motherindex2=fTR->genInfoMo2[i];
      if(motherindex>=0) IsParticleFromB(motherindex);
      if(motherindex2>=0) IsParticleFromB(motherindex2);
      if(nEvent.HasSoftLepton) {
	nEvent.SoftLeptonPt=fTR->genInfoPt[i];
	break;
      }
    }//end of isLepton
  }//end of gen particle loop
}

float GetCoreResolutionScalingFactor(const float jeta) {
  float ajeta=abs(jeta);
  if(ajeta<=1.1) return 1.07;
  if(ajeta<1.7) return 1.10;
  if(ajeta<2.3) return 1.07;
  if(ajeta<5.0) return 1.18;
  return 1.0; // not defined out here - don't smear
}

int JZBAnalysis::FindGenJetIndex(const float jpt, const float jeta, const float jphi) {
  int matchedindex=-1;
  float mindr=999.99;
  for(int ijet=0;ijet<fTR->NGenJets;ijet++) {
    double dr=sqrt( (jeta-fTR->GenJetEta[ijet]) * (jeta-fTR->GenJetEta[ijet]) + (jphi-fTR->GenJetPhi[ijet])*(jphi-fTR->GenJetPhi[ijet]));
    if(dr>0.3) continue;
    
    if(dr>mindr) continue;
    
    mindr=dr;
    matchedindex=ijet;
  }
  
  return matchedindex;
}

float JZBAnalysis::smearedJetPt(const float jpt, const float jeta, const float jphi) {
  //jet resolution oversmearing as described here: 
  //https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
  if(fDataType_!= "mc") return jpt; // if we're dealing with data we don't smear!

  int genIndex = FindGenJetIndex(jpt,jeta,jphi);
  if(genIndex<0) return -jpt; // no associated gen jet
  float genPt = fTR->GenJetPt[genIndex];
  float c = GetCoreResolutionScalingFactor(jeta);
  return max((float)0.,genPt+c*(jpt-genPt));
}

bool JZBAnalysis::IsPUJet(const float jpt, const float jeta, const float jphi) {
  if(fDataType_!= "mc") return false; // if we're dealing with data we label all jets as not being PU related

  int genIndex = FindGenJetIndex(jpt,jeta,jphi);
  if(genIndex<0) return true; // no associated gen jet -> PU!
  return false;
}

  

void JZBAnalysis::addPath(std::vector<std::string>& paths, std::string base,
                          unsigned int start, unsigned int end) {

  for ( unsigned int i=start; i<=end; ++i ) {
    ostringstream path;
    path << base << "_v" << i;
    paths.push_back(path.str());
  }
  std::cout << "Added " << base << " (v" << start;
  if ( start!=end) std::cout << "-v" << end;
  std::cout << ")" << std::endl;

}

int split(string s, string ch, string temp[]) { 
  //splits a string into an array of strings using ch as delimiter
  int cnt=0, ndx;
  while (s!="") {  
    ndx=s.find(ch);
    if (ndx<0) {temp[cnt]=s; s=""; }
    else if (ndx==0) {temp[cnt]=""; s=s.substr(1); }
    else {temp[cnt]=s.substr(0,ndx);s=s.substr(ndx+1);  }
    cnt++;
  }
  return cnt;
}

int JZBAnalysis::ExtractFileNumber(string pathname) {
  if(pathname=="") return -1;
  size_t found=pathname.find_last_of("/\\");
  string file = pathname.substr(found+1);
  string Parts[20];
  for(int i=0;i<20;i++) Parts[i]="";
  split(file, "_", Parts);
  if(Parts[3]=="") return -1;
  return atoi(Parts[3].c_str()); // 1 is NTUpleProducer, 2 is CMSSW version, 3 is file number (!), 4 is retry number, 5 is "unique" identifier.root
}

bool JZBAnalysis::IsThisDY(vector<string> fileList) {
  bool isDY=false;
  for(int ifile=0;ifile<fileList.size();ifile++) {
    if((int)fileList[ifile].find("/DY")>-1) isDY=true;
  }
  if(isDY) return true;
  else return false;
}

JZBAnalysis::JZBAnalysis(TreeReader *tr, std::string dataType, bool fullCleaning, bool isModelScan, bool makeSmall, bool doGenInfo, vector<string> fileList) :
  UserAnalysisBase(tr,dataType!="mc"), fDataType_(dataType), fFullCleaning_(fullCleaning) , fisModelScan(isModelScan) , fmakeSmall(makeSmall), fdoGenInfo(doGenInfo) {
    
  if(fileList.size()==1) fFile=ExtractFileNumber(fileList[0]);
  else fFile=-1;
  fIsDY=IsThisDY(fileList);
  // Define trigger paths to check
  addPath(elTriggerPaths,"HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL",10,30);

  addPath(muTriggerPaths,"HLT_Mu13_Mu8",10,30);
  addPath(muTriggerPaths,"HLT_Mu17_Mu8",10,30);
  addPath(muTriggerPaths,"HLT_Mu17_TkMu8",9,20);
  
  addPath(meTriggerPaths,"HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL", 4, 15);
  addPath(emTriggerPaths,"HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL", 4, 15);

  addPath(metTriggerPaths,"HLT_PFMET150", 2, 10);
  
  //FIXME: to be filled in with info. from Kostas
  addPath(htTriggerPaths,"", 0, 0);

  addPath(singleElTriggerPaths, "HLT_Ele27_WP80", 0, 20);
  addPath(singleMuTriggerPaths, "HLT_IsoMu24_eta2p1", 0, 20);

  
  L5corr_bJ = new TF1("L5corr_bJ","[0]+log10((x-[3])/[4])*([1]+[2]*log10((x-[3])/[4]))",10,3650);
  L5corr_bJ->SetParameters(1.05724,-0.031178,0.009538,6.80258,0.851021); //values from /shome/buchmann/material/JEStxtfiles/GR_R_50_V9_L5Flavor_bJ_AK5PFchs.txt
  L5corr_cJ = new TF1("L5corr_cJ","[0]+log10((x-[3])/[4])*([1]+[2]*log10((x-[3])/[4]))",10,3650);
  L5corr_cJ->SetParameters(0.854123,0.086427,-0.010176,5.72599,0.93306); //values from /shome/buchmann/material/JEStxtfiles/GR_R_50_V9_L5Flavor_cJ_AK5PFchs.txt
  L5corr_qJ = new TF1("L5corr_qJ","[0]+log10((x-[3])/[4])*([1]+[2]*log10((x-[3])/[4]))",10,3650);
  L5corr_qJ->SetParameters(0.691709,0.176237,-0.026105,4.94774,1.00526); //values from /shome/buchmann/material/JEStxtfiles/GR_R_50_V9_L5Flavor_qJ_AK5PFchs.txt
  L5corr_gJ = new TF1("L5corr_gJ","[0]+log10((x-[3])/[4])*([1]+[2]*log10((x-[3])/[4]))",10,3650);
  L5corr_gJ->SetParameters(1.08381,-0.044501,0.006947,3.00346,0.974058); //values from /shome/buchmann/material/JEStxtfiles/GR_R_50_V9_L5Flavor_gJ_AK5PFchs.txt
}

//________________________________________________________________________________________
const bool JZBAnalysis::passFilters( int& bits ) {

  // Check event filters
  bits = 0;
  if ( !fTR->HBHENoiseFilterResult ) bits |= 1;
  if ( !fTR->hcalLaserEventFilter )  bits |= (1<<2);
  if ( !fTR->EcalDeadCellTriggerPrimitiveFilter ) bits |= (1<<3);
  if ( !fTR->trackingFailureFilter ) bits |= (1<<4);
  if ( !fTR->eeBadScFilter )         bits |= (1<<5);
  if ( !fTR->CSCTightHaloID )        bits |= (1<<6);

  return (bits==0);

}

//________________________________________________________________________________________
JZBAnalysis::~JZBAnalysis(){}

//________________________________________________________________________________________
void JZBAnalysis::Begin(TFile *f){

  
  CSVT_CorrectionFile = new TFile("/shome/buchmann/material/Corrections/CSVT_Complete.root");
  CSVT_EfficiencyCorrection = (TH2F*)CSVT_CorrectionFile->Get("EfficiencyCorrection");
  CSVT_EfficiencyCorrectionUncert = (TH2F*)CSVT_CorrectionFile->Get("EfficiencyCorrectionUncertainty");
  CSVT_MisTagCorrection = (TH2F*)CSVT_CorrectionFile->Get("MisTag");
  CSVT_MisTagCorrectionUncert = (TH2F*)CSVT_CorrectionFile->Get("MisTagUncertainty");

  CSVM_CorrectionFile = new TFile("/shome/buchmann/material/Corrections/CSVM_Complete.root");
  CSVM_EfficiencyCorrection = (TH2F*)CSVM_CorrectionFile->Get("EfficiencyCorrection");
  CSVM_EfficiencyCorrectionUncert = (TH2F*)CSVM_CorrectionFile->Get("EfficiencyCorrectionUncertainty");
  CSVM_MisTagCorrection = (TH2F*)CSVM_CorrectionFile->Get("MisTag");
  CSVM_MisTagCorrectionUncert = (TH2F*)CSVM_CorrectionFile->Get("MisTagUncertainty");

  CSVL_CorrectionFile = new TFile("/shome/buchmann/material/Corrections/CSVL_Complete.root");
  CSVL_EfficiencyCorrection = (TH2F*)CSVL_CorrectionFile->Get("EfficiencyCorrection");
  CSVL_EfficiencyCorrectionUncert = (TH2F*)CSVL_CorrectionFile->Get("EfficiencyCorrectionUncertainty");
  CSVL_MisTagCorrection = (TH2F*)CSVL_CorrectionFile->Get("MisTag");
  CSVL_MisTagCorrectionUncert = (TH2F*)CSVL_CorrectionFile->Get("MisTagUncertainty");
  
  f->cd();

  rand_ = new TRandom();
  r = new TRandom3();
  
  TH1::AddDirectory(kFALSE);
  myInfo = new TTree("info","info/S");
  TString *user = new TString();
  TString *timestamp = new TString();
  TString *jzbversion = new TString();
  TString *cmsdir = new TString();
  TString *jzbinfo = new TString();
  myInfo->Branch("user",&user,16000,0);
  myInfo->Branch("timestamp",&timestamp,16000,0);
  myInfo->Branch("version",&jzbversion,16000,0);
  myInfo->Branch("cmsdir",&cmsdir,16000,0);
  myInfo->Branch("jzbinfo",&jzbinfo,16000,0);
  char usertext[255];
  *jzbversion=sjzbversion;
  char scmsdir[1000];
  getcwd(scmsdir,1000);
  *cmsdir=scmsdir;
  *jzbinfo=sjzbinfo;
  *user=getenv("USER");
  time_t rawtime;
  struct tm * timeinfo;
  time (&rawtime );
  *timestamp=ctime(&rawtime);
  
  f->cd();
  myInfo->Fill();
  myInfo->Write();

  weight_histo = new TH1F("weight_histo","weight_histo",1,0,2);
  
  myTree = new TTree("events","events");

  myTree->Branch("NTupleNumber",&fFile,"NTupleNumber/I");
  myTree->Branch("is_data",&nEvent.is_data,"is_data/O");
  myTree->Branch("mll",&nEvent.mll,"mll/F");
  myTree->Branch("minc",&nEvent.minc,"minc/F");
  myTree->Branch("pt",&nEvent.pt,"pt/F");
  myTree->Branch("phi",&nEvent.phi,"phi/F");
  myTree->Branch("eta",&nEvent.eta,"eta/F");
  myTree->Branch("E",&nEvent.E,"E/F");
  myTree->Branch("pt1",&nEvent.pt1,"pt1/F");
  myTree->Branch("pt2",&nEvent.pt2,"pt2/F");
  myTree->Branch("l1l2dR",&nEvent.l1l2dR,"l1l2dR/F");
  myTree->Branch("iso1",&nEvent.iso1,"iso1/F");
  myTree->Branch("iso2",&nEvent.iso2,"iso2/F");
  myTree->Branch("isConv1",&nEvent.isConv1,"isConv1/O");
  myTree->Branch("isConv2",&nEvent.isConv2,"isConv2/O");
  myTree->Branch("softMuon",&nEvent.softMuon,"softMuon/O");
  myTree->Branch("softMuonMC",&nEvent.softMuonMC,"softMuonMC/O");

  myTree->Branch("NgenZs",&nEvent.NgenZs,"NgenZs/I");
  myTree->Branch("NgenLeps",&nEvent.NgenLeps,"NgenLeps/I");
  
  myTree->Branch("genPt1",&nEvent.genPt1,"genPt1/F");
  myTree->Branch("genPt2",&nEvent.genPt2,"genPt2/F");
  myTree->Branch("genDRll",&nEvent.genDRll,"genDRll/F");
  myTree->Branch("genPhi1",&nEvent.genPhi1,"genPhi1/F");
  myTree->Branch("genPhi2",&nEvent.genPhi2,"genPhi2/F");
  myTree->Branch("genEta1",&nEvent.genEta1,"genEta1/F");
  myTree->Branch("genEta2",&nEvent.genEta2,"genEta2/F");
  myTree->Branch("genId1",&nEvent.genId1,"genId1/I");
  myTree->Branch("genId2",&nEvent.genId2,"genId2/I");
  myTree->Branch("genMID1gen",&nEvent.genMID1gen,"genMID1gen/I");
  myTree->Branch("genMID2gen",&nEvent.genMID2gen,"genMID2gen/I");
  myTree->Branch("genMID1",&nEvent.genMID1,"genMID1/I");
  myTree->Branch("genMID2",&nEvent.genMID2,"genMID2/I");
  myTree->Branch("genGMID1gen",&nEvent.genGMID1gen,"genGMID1gen/I");
  myTree->Branch("genGMID2gen",&nEvent.genGMID2gen,"genGMID2gen/I");
  myTree->Branch("genGMID1",&nEvent.genGMID1,"genGMID1/I");
  myTree->Branch("genGMID2",&nEvent.genGMID2,"genGMID2/I");
  myTree->Branch("genMET",&nEvent.genMET,"genMET/F");
  myTree->Branch("genZPt",&nEvent.genZPt,"genZPt/F");
  myTree->Branch("genMll",&nEvent.genMll,"genMll/F");
  myTree->Branch("genRecoil",&nEvent.genRecoil,"genRecoil/F");
  myTree->Branch("genJZB",&nEvent.genJZB,"genJZB/F");
  myTree->Branch("genNjets",&nEvent.genNjets,"genNjets/I");
  myTree->Branch("genNleptons",&nEvent.genNleptons,"genNleptons/I");
  myTree->Branch("genPt1Sel",&nEvent.genPt1Sel,"genPt1Sel/F");
  myTree->Branch("genPt2Sel",&nEvent.genPt2Sel,"genPt2Sel/F");
  myTree->Branch("genEta1Sel",&nEvent.genEta1Sel,"genEta1Sel/F");
  myTree->Branch("genEta2Sel",&nEvent.genEta2Sel,"genEta2Sel/F");
  myTree->Branch("genId1Sel",&nEvent.genId1Sel,"genId1Sel/I");
  myTree->Branch("genId2Sel",&nEvent.genId2Sel,"genId2Sel/I");
  myTree->Branch("genZPtSel",&nEvent.genZPtSel,"genZPtSel/F");
  myTree->Branch("genMllSel",&nEvent.genMllSel,"genMllSel/F");
  myTree->Branch("genRecoilSel",&nEvent.genRecoilSel,"genRecoilSel/F");
  myTree->Branch("genJZBSel",&nEvent.genJZBSel,"genJZBSel/F");
  myTree->Branch("eta1",&nEvent.eta1,"eta1/F");
  myTree->Branch("eta2",&nEvent.eta2,"eta2/F");
  myTree->Branch("phi1",&nEvent.phi1,"phi1/F");
  myTree->Branch("phi2",&nEvent.phi2,"phi2/F");
  myTree->Branch("dphiZpfMet",&nEvent.dphiZpfMet,"dphiZpfMet/F");
  myTree->Branch("dphiZs1",&nEvent.dphiZs1,"dphiZs1/F");
  myTree->Branch("dphiZs2",&nEvent.dphiZs2,"dphiZs2/F");
  myTree->Branch("dphiMet1",&nEvent.dphiMet1,"dphiMet1/F");
  myTree->Branch("dphiMet2",&nEvent.dphiMet2,"dphiMet2/F");
  myTree->Branch("dphitcMet1",&nEvent.dphitcMet1,"dphitcMet1/F");
  myTree->Branch("dphitcMet2",&nEvent.dphitcMet2,"dphitcMet2/F");
  myTree->Branch("dphipft1Met1",&nEvent.dphipft1Met1,"dphipft1Met1/F");
  myTree->Branch("dphipft1Met2",&nEvent.dphipft1Met2,"dphipft1Met2/F");
  myTree->Branch("dphipfRecoilMet1",&nEvent.dphipfRecoilMet1,"dphipfRecoilMet1/F");
  myTree->Branch("dphipfRecoilMet2",&nEvent.dphipfRecoilMet2,"dphipfRecoilMet2/F");
  myTree->Branch("dphi",&nEvent.dphi,"dphi/F");
  myTree->Branch("ElCInfoIsGsfCtfCons",&nEvent.ElCInfoIsGsfCtfCons,"ElCInfoIsGsfCtfCons/O");
  myTree->Branch("ElCInfoIsGsfScPixCons",&nEvent.ElCInfoIsGsfScPixCons,"ElCInfoIsGsfScPixCons/O");
  myTree->Branch("ElCInfoIsGsfCtfScPixCons",&nEvent.ElCInfoIsGsfCtfScPixCons,"ElCInfoIsGsfCtfScPixCons/O");

  myTree->Branch("id1",&nEvent.id1,"id1/I");
  myTree->Branch("id2",&nEvent.id2,"id2/I");
  myTree->Branch("ch1",&nEvent.ch1,"ch1/I");
  myTree->Branch("ch2",&nEvent.ch2,"ch2/I");
  myTree->Branch("chid1",&nEvent.chid1,"chid1/I");
  myTree->Branch("chid2",&nEvent.chid2,"chid2/I");
  myTree->Branch("process",&nEvent.process,"process/I");

  myTree->Branch("leptonNum",&nEvent.leptonNum,"leptonNum/I");
  myTree->Branch("leptonPt",nEvent.leptonPt,"leptonPt[leptonNum]/F");
  myTree->Branch("leptonEta",nEvent.leptonEta,"leptonEta[leptonNum]/F");
  myTree->Branch("leptonPhi",nEvent.leptonPhi,"leptonPhi[leptonNum]/F");
  myTree->Branch("leptonId",nEvent.leptonId,"leptonId[leptonNum]/I");
  myTree->Branch("leptonCharge",nEvent.leptonCharge,"leptonCharge[leptonNum]/I");

  myTree->Branch("leptonPairNum",&nEvent.leptonPairNum,"leptonPairNum/I");
  myTree->Branch("leptonPairMass",nEvent.leptonPairMass,"leptonPairMass[leptonPairNum]/F");
  myTree->Branch("leptonPairDphi",nEvent.leptonPairDphi,"leptonPairDphi[leptonPairNum]/F");
  myTree->Branch("leptonPairId",nEvent.leptonPairId,"leptonPairId[leptonPairNum]/I");

  myTree->Branch("met",nEvent.met,Form("met[%d]/F",int(mettype_max)));
  myTree->Branch("pfMET",&nEvent.pfMET,"pfMET/F");
  myTree->Branch("t1pfMET",&nEvent.t1pfMET,"t1pfMET/F");
  myTree->Branch("metPhi",nEvent.metPhi,Form("metPhi[%d]/F",int(mettype_max)));
  myTree->Branch("dphiMetLep",nEvent.dphiMetLep,Form("dphiMetLep[%d]/F",int(mettype_max)));
  myTree->Branch("dphiMetJet",nEvent.dphiMetJet,Form("dphiMetJet[%d]/F",int(mettype_max)));
  myTree->Branch("dphiMetSumJetPt",nEvent.dphiMetSumJetPt,Form("dphiMetSumJetPt[%d]/F",int(mettype_max)));
  myTree->Branch("metPerp",nEvent.metPerp,Form("metPerp[%d]/F",int(mettype_max)));
  myTree->Branch("metPar",nEvent.metPar,Form("metPar[%d]/F",int(mettype_max)));


  myTree->Branch("eventNum",&nEvent.eventNum,"eventNum/l");
  myTree->Branch("runNum",&nEvent.runNum,"runNum/I");
  myTree->Branch("lumi",&nEvent.lumi,"lumi/I");
  myTree->Branch("goodVtx",&nEvent.goodVtx,"goodVtx/I");
  myTree->Branch("numVtx",&nEvent.numVtx,"numVtx/I");
  myTree->Branch("badJet",&nEvent.badJet,"badJet/I");
  myTree->Branch("totEvents",&nEvent.totEvents,"totEvents/F");

  myTree->Branch("pfJetNum",&nEvent.pfJetNum,"pfJetNum/I");
  myTree->Branch("pfJetPt",nEvent.pfJetPt,"pfJetPt[pfJetNum]/F");
  myTree->Branch("pfJetEta",nEvent.pfJetEta,"pfJetEta[pfJetNum]/F");
  myTree->Branch("pfJetPhi",nEvent.pfJetPhi,"pfJetPhi[pfJetNum]/F");
  myTree->Branch("pfJetID",nEvent.pfJetID,"pfJetID[pfJetNum]/O");
  myTree->Branch("pfJetScale",nEvent.pfJetScale,"pfJetScale[pfJetNum]/F");
  myTree->Branch("pfJetScaleUnc",nEvent.pfJetScaleUnc,"pfJetScaleUnc[pfJetNum]/F");
  myTree->Branch("pfJetDphiMet",nEvent.pfJetDphiMet,"pfJetDphiMet[pfJetNum]/F");
  myTree->Branch("pfJetDphiZ",nEvent.pfJetDphiZ,"pfJetDphiZ[pfJetNum]/F");
  myTree->Branch("pfBJetDphiZ",nEvent.pfBJetDphiZ,"pfBJetDphiZ[pfJetNum]/F");
  myTree->Branch("pfHT",&nEvent.pfHT,"pfHT/F");
  myTree->Branch("pfGoodHT",&nEvent.pfGoodHT,"pfGoodHT/F");
  myTree->Branch("pfTightHT",&nEvent.pfTightHT,"pfTightHT/F");
  myTree->Branch("CorrectionRatio", nEvent.CorrectionRatio,"CorrectionRatio[pfJetNum]/F");

  myTree->Branch("metUncertainty",&nEvent.metUncertainty,"metUncertainty/F");
  myTree->Branch("type1metUncertainty",&nEvent.type1metUncertainty,"type1metUncertainty/F");

  myTree->Branch("pfJetGoodNum40CHS",&nEvent.pfJetGoodNum40CHS,"pfJetGoodNum40CHS/I");
  myTree->Branch("pfJetGoodNum30Fwd",&nEvent.pfJetGoodNum30Fwd,"pfJetGoodNum30Fwd/I");
  myTree->Branch("pfJetGoodNum40Fwd",&nEvent.pfJetGoodNum40Fwd,"pfJetGoodNum40Fwd/I");
  myTree->Branch("pfJetGoodNum50Fwd",&nEvent.pfJetGoodNum50Fwd,"pfJetGoodNum50Fwd/I");
  myTree->Branch("pfJetGoodNum60Fwd",&nEvent.pfJetGoodNum60Fwd,"pfJetGoodNum60Fwd/I");
  myTree->Branch("pfJetGoodNum30",&nEvent.pfJetGoodNum30,"pfJetGoodNum30/I");
  myTree->Branch("pfJetGoodNum40",&nEvent.pfJetGoodNum40,"pfJetGoodNum40/I");
  myTree->Branch("pfJetGoodNum50",&nEvent.pfJetGoodNum50,"pfJetGoodNum50/I");
  myTree->Branch("pfJetGoodNum60",&nEvent.pfJetGoodNum60,"pfJetGoodNum60/I");
  myTree->Branch("pfJetGoodNumBtag30",&nEvent.pfJetGoodNumBtag30,"pfJetGoodNumBtag30/I");
  myTree->Branch("pfJetGoodNumBtag40",&nEvent.pfJetGoodNumBtag40,"pfJetGoodNumBtag40/I");
  myTree->Branch("pfJetGoodNumBtag30_Tight",&nEvent.pfJetGoodNumBtag30_Tight,"pfJetGoodNumBtag30_Tight/I");
  myTree->Branch("pfJetGoodNumBtag40_Tight",&nEvent.pfJetGoodNumBtag40_Tight,"pfJetGoodNumBtag40_Tight/I");
  myTree->Branch("pfJetGoodNumBtag30_Loose",&nEvent.pfJetGoodNumBtag30_Loose,"pfJetGoodNumBtag30_Loose/I");
  myTree->Branch("pfJetGoodNumBtag40_Loose",&nEvent.pfJetGoodNumBtag40_Loose,"pfJetGoodNumBtag40_Loose/I");
  myTree->Branch("pfJetGoodNumIDBtag",&nEvent.pfJetGoodNumIDBtag,"pfJetGoodNumIDBtag/I");
  myTree->Branch("pfJetGoodNumID",&nEvent.pfJetGoodNumID,"pfJetGoodNumID/I");
  myTree->Branch("pfJetGoodNump1sigma",&nEvent.pfJetGoodNump1sigma,"pfJetGoodNump1sigma/I");
  myTree->Branch("pfJetGoodNumn1sigma",&nEvent.pfJetGoodNumn1sigma,"pfJetGoodNumn1sigma/I");
  myTree->Branch("pfJetGoodNum40p1sigma",&nEvent.pfJetGoodNum40p1sigma,"pfJetGoodNum40p1sigma/I");
  myTree->Branch("pfJetGoodNum40n1sigma",&nEvent.pfJetGoodNum40n1sigma,"pfJetGoodNum40n1sigma/I");
  myTree->Branch("pfJetGoodNum50p1sigma",&nEvent.pfJetGoodNum50p1sigma,"pfJetGoodNum50p1sigma/I");
  myTree->Branch("pfJetGoodNum50n1sigma",&nEvent.pfJetGoodNum50n1sigma,"pfJetGoodNum50n1sigma/I");

  myTree->Branch("pfJetGoodPt", nEvent.pfJetGoodPt,"pfJetGoodPt[pfJetGoodNum40]/F");
  myTree->Branch("pfJetGoodTracks", nEvent.pfJetGoodTracks,"pfJetGoodTracks[pfJetGoodNum40]/I");
  myTree->Branch("pfJetGoodTracksN", nEvent.pfJetGoodTracksN,"pfJetGoodTracksN[pfJetGoodNum40]/I");
  myTree->Branch("pfJetGoodEta",nEvent.pfJetGoodEta,"pfJetGoodEta[pfJetGoodNum40]/F");
  myTree->Branch("pfJetGoodPhi",nEvent.pfJetGoodPhi,"pfJetGoodPhi[pfJetGoodNum40]/F");
  myTree->Branch("pfJetGoodID", nEvent.pfJetGoodID,"pfJetGoodID[pfJetGoodNum40]/F");
  myTree->Branch("pfJetGoodE", nEvent.pfJetGoodE,"pfJetGoodE[pfJetGoodNum40]/F");
  myTree->Branch("pfJetGoodM", nEvent.pfJetGoodM,"pfJetGoodM[pfJetGoodNum40]/F");
  myTree->Branch("pfJetGoodMl", nEvent.pfJetGoodMl,"pfJetGoodMl[pfJetGoodNum40]/F");
  myTree->Branch("pfJetGoodPtl", nEvent.pfJetGoodPtl,"pfJetGoodPtl[pfJetGoodNum40]/F");

  myTree->Branch("pfJetGoodPtBtag", nEvent.pfJetGoodPtBtag,"pfJetGoodPtBag[pfJetGoodNumBtag40]/F");
  myTree->Branch("pfJetGoodTracksBtag", nEvent.pfJetGoodTracksBtag,"pfJetGoodTracksBtag[pfJetGoodNumBtag40]/I");
  myTree->Branch("pfJetGoodTracksNBtag", nEvent.pfJetGoodTracksNBtag,"pfJetGoodTracksNBtag[pfJetGoodNumBtag40]/I");
  myTree->Branch("pfJetGoodEtaBtag",nEvent.pfJetGoodEtaBtag,"pfJetGoodEtaBtag[pfJetGoodNumBtag40]/F");
  myTree->Branch("pfJetGoodPhiBtag",nEvent.pfJetGoodPhiBtag,"pfJetGoodPhiBtag[pfJetGoodNumBtag40]/F");
  myTree->Branch("pfJetGoodIDBtag", nEvent.pfJetGoodIDBtag,"pfJetGoodID[pfJetGoodNumBtag40]/O");
  myTree->Branch("pfJetGoodEBtag", nEvent.pfJetGoodEBtag,"pfJetGoodEBtag[pfJetGoodNumBtag40]/F");
  myTree->Branch("pfJetGoodMBtag", nEvent.pfJetGoodMBtag,"pfJetGoodMBtag[pfJetGoodNumBtag40]/F");
  myTree->Branch("pfJetGoodMlBtag", nEvent.pfJetGoodMlBtag,"pfJetGoodMlBtag[pfJetGoodNumBtag40]/F");
  myTree->Branch("pfJetGoodPtlBtag", nEvent.pfJetGoodPtlBtag,"pfJetGoodPtlBtag[pfJetGoodNumBtag40]/F");
  myTree->Branch("bTagProbCSVBP", nEvent.bTagProbCSVBP,"bTagProbCSVBP[pfJetGoodNum40]/F");
  myTree->Branch("bTagProbCSVMVA", nEvent.bTagProbCSVMVA,"bTagProbCSVMVA[pfJetGoodNum40]/F");


  myTree->Branch("jzb",nEvent.jzb,Form("jzb[%d]/F",int(jzbtype_max)));
  myTree->Branch("sjzb",nEvent.sjzb,Form("sjzb[%d]/F",int(jzbtype_max)));
  myTree->Branch("dphi_sumJetVSZ",nEvent.dphi_sumJetVSZ,Form("dphi_sumJetVSZ[%d]/F",int(jzbtype_max)));
  myTree->Branch("sumJetPt",nEvent.sumJetPt,Form("sumJetPt[%d]/F",int(jzbtype_max)));

  myTree->Branch("weight", &nEvent.weight,"weight/F");
  myTree->Branch("PUweight",&nEvent.PUweight,"PUweight/F");
  myTree->Branch("PUweightUP",&nEvent.PUweightUP,"PUweightUP/F");
  myTree->Branch("PUweightDOWN",&nEvent.PUweightDOWN,"PUweightDOWN/F");
  myTree->Branch("Efficiencyweightonly",&nEvent.Efficiencyweightonly,"Efficiencyweightonly/F");
  myTree->Branch("weightEffDown",&nEvent.weightEffDown,"weightEffDown/F");
  myTree->Branch("weightEffUp",&nEvent.weightEffUp,"weightEffUp/F");

  
  myTree->Branch("ZbCHS3010_LeadingJetIsPu",&nEvent.ZbCHS3010_LeadingJetIsPu,"ZbCHS3010_LeadingJetIsPu/O");
  myTree->Branch("ZbCHS3010_SubLeadingJetIsPu",&nEvent.ZbCHS3010_SubLeadingJetIsPu,"ZbCHS3010_SubLeadingJetIsPu/O");
  
  myTree->Branch("ZbCHS1010_LeadingJetIsPu",&nEvent.ZbCHS1010_LeadingJetIsPu,"ZbCHS1010_LeadingJetIsPu/O");
  myTree->Branch("ZbCHS1010_SubLeadingJetIsPu",&nEvent.ZbCHS1010_SubLeadingJetIsPu,"ZbCHS1010_SubLeadingJetIsPu/O");

  myTree->Branch("ZbCHS3010_BTagWgtT",&nEvent.ZbCHS3010_BTagWgtT,"ZbCHS3010_BTagWgtT/F");
  myTree->Branch("ZbCHS1010_BTagWgtT",&nEvent.ZbCHS1010_BTagWgtT,"ZbCHS1010_BTagWgtT/F");
  
  myTree->Branch("ZbCHS3010_BTagWgtTUp",&nEvent.ZbCHS3010_BTagWgtTUp,"ZbCHS3010_BTagWgtTUp/F");
  myTree->Branch("ZbCHS1010_BTagWgtTUp",&nEvent.ZbCHS1010_BTagWgtTUp,"ZbCHS1010_BTagWgtTUp/F");
  
  myTree->Branch("ZbCHS3010_BTagWgtTDown",&nEvent.ZbCHS3010_BTagWgtTDown,"ZbCHS3010_BTagWgtTDown/F");
  myTree->Branch("ZbCHS1010_BTagWgtTDown",&nEvent.ZbCHS1010_BTagWgtTDown,"ZbCHS1010_BTagWgtTDown/F");
  
  myTree->Branch("ZbCHS3010_BTagWgtM",&nEvent.ZbCHS3010_BTagWgtM,"ZbCHS3010_BTagWgtM/F");
  myTree->Branch("ZbCHS1010_BTagWgtM",&nEvent.ZbCHS1010_BTagWgtM,"ZbCHS1010_BTagWgtM/F");
  
  myTree->Branch("ZbCHS3010_BTagWgtMUp",&nEvent.ZbCHS3010_BTagWgtMUp,"ZbCHS3010_BTagWgtMUp/F");
  myTree->Branch("ZbCHS1010_BTagWgtMUp",&nEvent.ZbCHS1010_BTagWgtMUp,"ZbCHS1010_BTagWgtMUp/F");
  
  myTree->Branch("ZbCHS3010_BTagWgtM",&nEvent.ZbCHS3010_BTagWgtMDown,"ZbCHS3010_BTagWgtMDown/F");
  myTree->Branch("ZbCHS1010_BTagWgtM",&nEvent.ZbCHS1010_BTagWgtMDown,"ZbCHS1010_BTagWgtMDown/F");
  
  myTree->Branch("ZbCHS3010_BTagWgtL",&nEvent.ZbCHS3010_BTagWgtL,"ZbCHS3010_BTagWgtL/F");
  myTree->Branch("ZbCHS1010_BTagWgtL",&nEvent.ZbCHS1010_BTagWgtL,"ZbCHS1010_BTagWgtL/F");
  
  myTree->Branch("ZbCHS3010_BTagWgtLUp",&nEvent.ZbCHS3010_BTagWgtLUp,"ZbCHS3010_BTagWgtLUp/F");
  myTree->Branch("ZbCHS1010_BTagWgtLUp",&nEvent.ZbCHS1010_BTagWgtLUp,"ZbCHS1010_BTagWgtLUp/F");
  
  myTree->Branch("ZbCHS3010_BTagWgtL",&nEvent.ZbCHS3010_BTagWgtLDown,"ZbCHS3010_BTagWgtLDown/F");
  myTree->Branch("ZbCHS1010_BTagWgtL",&nEvent.ZbCHS1010_BTagWgtLDown,"ZbCHS1010_BTagWgtLDown/F");
  
  myTree->Branch("passed_triggers", &nEvent.passed_triggers,"passed_triggers/O");
  myTree->Branch("trigger_bit", &nEvent.trigger_bit,"trigger_bit/I");
  myTree->Branch("passed_filters", &nEvent.passed_filters,"passed_filters/O");
  myTree->Branch("filter_bit", &nEvent.filter_bit,"filter_bit/I");

  myTree->Branch("MassGlu",&nEvent.mGlu,"MassGlu/F");
  myTree->Branch("MassChi",&nEvent.mChi,"MassChi/F");
  myTree->Branch("MassLSP",&nEvent.mLSP,"MassLSP/F");
  myTree->Branch("xSMS",&nEvent.xSMS,"xSMS/F");
  myTree->Branch("xbarSMS",&nEvent.xbarSMS,"xbarSMS/F");
  myTree->Branch("MassGMSBGlu",&nEvent.mGMSBGlu,"MassGlu/F");
  myTree->Branch("MassGMSBChi",&nEvent.mGMSBChi,"MassChi/F");
  myTree->Branch("MassGMSBLSP",&nEvent.mGMSBLSP,"MassLSP/F");
  myTree->Branch("NPdfs",&nEvent.NPdfs,"NPdfs/I");
  myTree->Branch("pdfW",nEvent.pdfW,"pdfW[NPdfs]/F");
  myTree->Branch("pdfWsum",&nEvent.pdfWsum,"pdfWsum/F");
  myTree->Branch("EventFlavor",&nEvent.EventFlavor,"EventFlavor/I");
  myTree->Branch("EventZToTaus",&nEvent.EventZToTaus,"EventZToTaus/O");

  myTree->Branch("ZbCHS3010_pfJetGoodNum",&nEvent.ZbCHS3010_pfJetGoodNum,"ZbCHS3010_pfJetGoodNum/I");
  myTree->Branch("ZbCHS3010_pfJetGoodNumBtag",&nEvent.ZbCHS3010_pfJetGoodNumBtag,"ZbCHS3010_pfJetGoodNumBtag/I");
  myTree->Branch("ZbCHS3010_bTagProbCSVBP",nEvent.ZbCHS3010_bTagProbCSVBP,"ZbCHS3010_bTagProbCSVBP[ZbCHS3010_pfJetGoodNum]/F");
  myTree->Branch("ZbCHS3010_pfJetGoodEta",nEvent.ZbCHS3010_pfJetGoodEta,"ZbCHS3010_pfJetGoodEta[ZbCHS3010_pfJetGoodNum]/F");
  myTree->Branch("ZbCHS3010_pfJetDphiZ",nEvent.ZbCHS3010_pfJetDphiZ,"ZbCHS3010_pfJetDphiZ[ZbCHS3010_pfJetGoodNum]/F");
  myTree->Branch("ZbCHS3010_pfJetGoodPt",nEvent.ZbCHS3010_pfJetGoodPt,"ZbCHS3010_pfJetGoodPt[ZbCHS3010_pfJetGoodNum]/F");
  myTree->Branch("ZbCHS3010_pfBJetDphiZ",nEvent.ZbCHS3010_pfBJetDphiZ,"ZbCHS3010_pfBJetDphiZ[ZbCHS3010_pfJetGoodNumBtag]/F");

  myTree->Branch("ZbCHS1010_pfJetGoodNum",&nEvent.ZbCHS1010_pfJetGoodNum,"ZbCHS1010_pfJetGoodNum/I");
  myTree->Branch("ZbCHS1010_pfJetGoodNumBtag",&nEvent.ZbCHS1010_pfJetGoodNumBtag,"ZbCHS1010_pfJetGoodNumBtag/I");
  myTree->Branch("ZbCHS1010_bTagProbCSVBP",nEvent.ZbCHS1010_bTagProbCSVBP,"ZbCHS1010_bTagProbCSVBP[ZbCHS1010_pfJetGoodNum]/F");
  myTree->Branch("ZbCHS1010_pfJetGoodEta",nEvent.ZbCHS1010_pfJetGoodEta,"ZbCHS1010_pfJetGoodEta[ZbCHS1010_pfJetGoodNum]/F");
  myTree->Branch("ZbCHS1010_pfJetDphiZ",nEvent.ZbCHS1010_pfJetDphiZ,"ZbCHS1010_pfJetDphiZ[ZbCHS1010_pfJetGoodNum]/F");
  myTree->Branch("ZbCHS1010_pfJetGoodPt",nEvent.ZbCHS1010_pfJetGoodPt,"ZbCHS1010_pfJetGoodPt[ZbCHS1010_pfJetGoodNum]/F");
  myTree->Branch("ZbCHS1010_pfBJetDphiZ",nEvent.ZbCHS1010_pfBJetDphiZ,"ZbCHS1010_pfBJetDphiZ[ZbCHS1010_pfJetGoodNumBtag]/F");
    
  myTree->Branch("MetFactor",&nEvent.fact,"MetFactor/F");
  myTree->Branch("tri_pt1",&nEvent.tri_pt1,"tri_pt1/F");
  myTree->Branch("tri_pt2",&nEvent.tri_pt2,"tri_pt2/F");
  myTree->Branch("tri_pt3",&nEvent.tri_pt3,"tri_pt3/F");
  myTree->Branch("tri_eta1",&nEvent.tri_eta1,"tri_eta1/F");
  myTree->Branch("tri_eta2",&nEvent.tri_eta2,"tri_eta2/F");
  myTree->Branch("tri_eta3",&nEvent.tri_eta3,"tri_eta3/F");
  myTree->Branch("tri_MatchFound",&nEvent.tri_MatchFound,"tri_MatchFound/O");
  myTree->Branch("tri_id1",&nEvent.tri_id1,"tri_id1/I");
  myTree->Branch("tri_id2",&nEvent.tri_id2,"tri_id2/I");
  myTree->Branch("tri_id3",&nEvent.tri_id3,"tri_id3/I");
  myTree->Branch("tri_ch1",&nEvent.tri_ch1,"tri_ch1/I");
  myTree->Branch("tri_ch2",&nEvent.tri_ch2,"tri_ch2/I");
  myTree->Branch("tri_ch3",&nEvent.tri_ch3,"tri_ch3/I");
  myTree->Branch("tri_mlll",&nEvent.tri_mlll,"tri_mlll/F");
  myTree->Branch("tri_mll",&nEvent.tri_mll,"tri_mll/F");
  myTree->Branch("tri_submll",&nEvent.tri_submll,"tri_submll/F");
  myTree->Branch("tri_mT",&nEvent.tri_mT,"tri_mT/F");
  myTree->Branch("tri_index1",&nEvent.tri_index1,"tri_index1/I");
  myTree->Branch("tri_index2",&nEvent.tri_index2,"tri_index2/I");
  myTree->Branch("tri_index3",&nEvent.tri_index3,"tri_index3/I");
  
  myTree->Branch("tri_badsubmll",&nEvent.tri_badsubmll,"tri_badsubmll/F");
  myTree->Branch("tri_badmll",&nEvent.tri_badmll,"tri_badmll/F");
  myTree->Branch("tri_badmT",&nEvent.tri_badmT,"tri_badmT/F");
  myTree->Branch("tri_badgenMID1",&nEvent.tri_badgenMID1,"tri_badgenMID1/I");
  myTree->Branch("tri_badgenMID2",&nEvent.tri_badgenMID2,"tri_badgenMID2/I");
  myTree->Branch("tri_badgenMID3",&nEvent.tri_badgenMID3,"tri_badgenMID3/I");
  myTree->Branch("tri_GoodZMatch",&nEvent.tri_GoodZMatch,"tri_GoodZMatch/O");
  myTree->Branch("tri_GoodWMatch",&nEvent.tri_GoodWMatch,"tri_GoodWMatch/O");
  
  myTree->Branch("tri_genMID1",&nEvent.tri_genMID1,"tri_genMID1/I");
  myTree->Branch("tri_genMID2",&nEvent.tri_genMID2,"tri_genMID2/I");
  myTree->Branch("tri_genMID3",&nEvent.tri_genMID3,"tri_genMID3/I");
  
  
  myTree->Branch("gentri_badsubmll",&nEvent.gentri_badsubmll,"gentri_badsubmll/F");
  myTree->Branch("gentri_badmll",&nEvent.gentri_badmll,"gentri_badmll/F");
  myTree->Branch("gentri_badmT",&nEvent.gentri_badmT,"gentri_badmT/F");
  myTree->Branch("gentri_badgenMID1",&nEvent.gentri_badgenMID1,"gentri_badgenMID1/I");
  myTree->Branch("gentri_badgenMID2",&nEvent.gentri_badgenMID2,"gentri_badgenMID2/I");
  myTree->Branch("gentri_badgenMID3",&nEvent.gentri_badgenMID3,"gentri_badgenMID3/I");
  myTree->Branch("gentri_GoodZMatch",&nEvent.gentri_GoodZMatch,"gentri_GoodZMatch/O");
  myTree->Branch("gentri_GoodWMatch",&nEvent.gentri_GoodWMatch,"gentri_GoodWMatch/O");

  myTree->Branch("tri_badindex1",&nEvent.tri_badindex1,"tri_badindex1/I");
  myTree->Branch("tri_badindex2",&nEvent.tri_badindex2,"tri_badindex2/I");
  myTree->Branch("tri_badindex3",&nEvent.tri_badindex3,"tri_badindex3/I");

  myTree->Branch("tri_badid1",&nEvent.tri_badid1,"tri_badid1/I");
  myTree->Branch("tri_badid2",&nEvent.tri_badid2,"tri_badid2/I");
  myTree->Branch("tri_badid3",&nEvent.tri_badid3,"tri_badid3/I");

  myTree->Branch("tri_dR12",&nEvent.tri_dR12,"tri_dR12/F");
  myTree->Branch("tri_dR13",&nEvent.tri_dR13,"tri_dR13/F");
  myTree->Branch("tri_dR23",&nEvent.tri_dR23,"tri_dR23/F");

  myTree->Branch("gentri_eta1",&nEvent.gentri_eta1,"gentri_eta1/F");
  myTree->Branch("gentri_eta",&nEvent.gentri_eta2,"gentri_eta2/F");
  myTree->Branch("gentri_eta3",&nEvent.gentri_eta3,"gentri_eta3/F");  
  myTree->Branch("gentri_pt1",&nEvent.gentri_pt1,"gentri_pt1/F");
  myTree->Branch("gentri_pt2",&nEvent.gentri_pt2,"gentri_pt2/F");
  myTree->Branch("gentri_pt3",&nEvent.gentri_pt3,"gentri_pt3/F");  
  myTree->Branch("gentri_id1",&nEvent.gentri_id1,"gentri_id1/F");
  myTree->Branch("gentri_id2",&nEvent.gentri_id2,"gentri_id2/F");
  myTree->Branch("gentri_id3",&nEvent.gentri_id3,"gentri_id3/F");
  myTree->Branch("gentri_ch1",&nEvent.gentri_ch1,"gentri_ch1/F");
  myTree->Branch("gentri_ch2",&nEvent.gentri_ch2,"gentri_ch2/F");
  myTree->Branch("gentri_ch3",&nEvent.gentri_ch3,"gentri_ch3/F");
  myTree->Branch("gentri_mlll",&nEvent.gentri_mlll,"gentri_mlll/F");
  myTree->Branch("gentri_mll",&nEvent.gentri_mll,"gentri_mll/F");
  myTree->Branch("gentri_submll",&nEvent.gentri_submll,"gentri_submll/F");
  myTree->Branch("gentri_mT",&nEvent.gentri_mT,"gentri_mT/F");
  
  myTree->Branch("gentri_genMID1",&nEvent.gentri_genMID1,"gentri_genMID1/I");
  myTree->Branch("gentri_genMID2",&nEvent.gentri_genMID2,"gentri_genMID2/I");
  myTree->Branch("gentri_genMID3",&nEvent.gentri_genMID3,"gentri_genMID3/I");
  
  myTree->Branch("gentri_dR12",&nEvent.gentri_dR12,"gentri_dR12/F");
  myTree->Branch("gentri_dR13",&nEvent.gentri_dR13,"gentri_dR13/F");
  myTree->Branch("gentri_dR23",&nEvent.gentri_dR23,"gentri_dR23/F");

  myTree->Branch("gentri_GoodZMatch",&nEvent.gentri_GoodZMatch,"gentri_GoodZMatch/O");
  myTree->Branch("gentri_GoodWMatch",&nEvent.gentri_GoodWMatch,"gentri_GoodWMatch/O");
  

  
  myTree->Branch("pgentri_eta1",&nEvent.pgentri_eta1,"pgentri_eta1/F");
  myTree->Branch("pgentri_eta",&nEvent.pgentri_eta2,"pgentri_eta2/F");
  myTree->Branch("pgentri_eta3",&nEvent.pgentri_eta3,"pgentri_eta3/F");  
  myTree->Branch("pgentri_pt1",&nEvent.pgentri_pt1,"pgentri_pt1/F");
  myTree->Branch("pgentri_pt2",&nEvent.pgentri_pt2,"pgentri_pt2/F");
  myTree->Branch("pgentri_pt3",&nEvent.pgentri_pt3,"pgentri_pt3/F");  
  myTree->Branch("pgentri_id1",&nEvent.pgentri_id1,"pgentri_id1/F");
  myTree->Branch("gentri_id2",&nEvent.pgentri_id2,"pgentri_id2/F");
  myTree->Branch("pgentri_id3",&nEvent.pgentri_id3,"pgentri_id3/F");
  myTree->Branch("pgentri_ch1",&nEvent.pgentri_ch1,"pgentri_ch1/F");
  myTree->Branch("pgentri_ch2",&nEvent.pgentri_ch2,"pgentri_ch2/F");
  myTree->Branch("pgentri_ch3",&nEvent.pgentri_ch3,"pgentri_ch3/F");
  myTree->Branch("pgentri_mlll",&nEvent.pgentri_mlll,"pgentri_mlll/F");
  myTree->Branch("pgentri_mll",&nEvent.pgentri_mll,"pgentri_mll/F");
  myTree->Branch("pgentri_submll",&nEvent.pgentri_submll,"pgentri_submll/F");
  myTree->Branch("pgentri_mT",&nEvent.pgentri_mT,"pgentri_mT/F");
  
  myTree->Branch("pgentri_dR12",&nEvent.pgentri_dR12,"pgentri_dR12/F");
  myTree->Branch("pgentri_dR13",&nEvent.pgentri_dR13,"pgentri_dR13/F");
  myTree->Branch("pgentri_dR23",&nEvent.pgentri_dR23,"pgentri_dR23/F");

  
  //generator information
  if(fdoGenInfo) {
	myTree->Branch("nZ",&nEvent.nZ,"nZ/I");
	myTree->Branch("SourceOfZ",&nEvent.SourceOfZ,"SourceOfZ[nZ]/I");
	myTree->Branch("DecayCode",&nEvent.DecayCode,"DecayCode/I");
	myTree->Branch("pureGeneratorJZB",&nEvent.pureGeneratorJZB,"pureGeneratorJZB/F");
	myTree->Branch("pureGeneratorZpt",&nEvent.pureGeneratorZpt,"pureGeneratorZpt/F");
	myTree->Branch("pureGeneratorZM",&nEvent.pureGeneratorZM,"pureGeneratorZM/F");
	myTree->Branch("pureGeneratorZeta",&nEvent.pureGeneratorZeta,"pureGeneratorZeta/F");
	myTree->Branch("pureGeneratorZphi",&nEvent.pureGeneratorZphi,"pureGeneratorZphi/F");
	myTree->Branch("pure2ndGeneratorJZB",&nEvent.pure2ndGeneratorJZB,"pure2ndGeneratorJZB/F");
	myTree->Branch("pure2ndGeneratorZpt",&nEvent.pure2ndGeneratorZpt,"pure2ndGeneratorZpt/F");
	myTree->Branch("LSPPromptnessLevel",&nEvent.LSPPromptnessLevel,"LSPPromptnessLevel[2]/I");
	myTree->Branch("ZPromptnessLevel",&nEvent.ZPromptnessLevel,"ZPromptnessLevel[2]/I");
	myTree->Branch("LSP1pt",&nEvent.LSP1pt,"LSP1pt/F");
	myTree->Branch("LSP2pt",&nEvent.LSP2pt,"LSP2pt/F");

	myTree->Branch("pureGeneratorMet",&nEvent.pureGeneratorMet,"pureGeneratorMet/F");
	myTree->Branch("pureGeneratorMetPhi",&nEvent.pureGeneratorMetPhi,"pureGeneratorMetPhi/F");
	myTree->Branch("pureGeneratorSumJetPt",&nEvent.pureGeneratorSumJetPt,"pureGeneratorSumJetPt/F");
	myTree->Branch("pureGeneratorSumJetEta",&nEvent.pureGeneratorSumJetEta,"pureGeneratorSumJetEta/F");
	myTree->Branch("pureGeneratorSumJetPhi",&nEvent.pureGeneratorSumJetPhi,"pureGeneratorSumJetPhi/F");

	myTree->Branch("LSP1Mo",&nEvent.LSP1Mo,"LSP1Mo/I");
	myTree->Branch("LSP2Mo",&nEvent.LSP2Mo,"LSP2Mo/I");
	myTree->Branch("LSP1Mopt",&nEvent.LSP1Mopt,"LSP1Mopt/F");
	myTree->Branch("LSP2Mopt",&nEvent.LSP2Mopt,"LSP2Mopt/F");

	myTree->Branch("nLSPs",&nEvent.nLSPs,"nLSPs/I");
	myTree->Branch("angleLSPLSP",&nEvent.angleLSPLSP,"angleLSPLSP/F");
	myTree->Branch("angleLSPLSP2d",&nEvent.angleLSPLSP2d,"angleLSPLSP2d/F");
	
	myTree->Branch("angleLSPZ2d",&nEvent.angleLSPZ2d,"angleLSPZ2d/F");
	myTree->Branch("angleChi2Z2d",&nEvent.angleChi2Z2d,"angleChi2Z2d/F");
	myTree->Branch("angleChi2Z",&nEvent.angleChi2Z,"angleChi2Z/F");
	myTree->Branch("angleLSPZ",&nEvent.angleLSPZ,"angleLSPZ/F");
	myTree->Branch("dphiSumLSPgenMET",&nEvent.dphiSumLSPgenMET,"dphiSumLSPgenMET/F");
	myTree->Branch("absvalSumLSP",&nEvent.absvalSumLSP,"absvalSumLSP/F");
	myTree->Branch("dphigenZgenMet",&nEvent.dphigenZgenMet,"dphigenZgenMet/F");
	myTree->Branch("SumLSPEta",&nEvent.SumLSPEta,"SumLSPEta/F");
	myTree->Branch("SumLSPPhi",&nEvent.SumLSPPhi,"SumLSPPhi/F");
  }

  myTree->Branch("realx",&nEvent.realx,"realx/F");
  myTree->Branch("imposedx",&nEvent.imposedx,"imposedx/F");
  
    //Z+b variables
  myTree->Branch("ZbCHS3010_alpha",&nEvent.ZbCHS3010_alpha,"ZbCHS3010_alpha/F");
  myTree->Branch("ZbCHS3010_alphaL5",&nEvent.ZbCHS3010_alphaL5,"ZbCHS3010_alphaL5/F");
  myTree->Branch("ZbCHS3010_L5corr",&nEvent.ZbCHS3010_L5corr,"ZbCHS3010_L5corr/F");
  myTree->Branch("ZbCHS1010_alpha",&nEvent.ZbCHS1010_alpha,"ZbCHS1010_alpha/F");
  myTree->Branch("ZbCHS3010_alphaUp",&nEvent.ZbCHS3010_alphaUp,"ZbCHS3010_alphaUp/F");
  myTree->Branch("ZbCHS1010_alphaUp",&nEvent.ZbCHS1010_alphaUp,"ZbCHS1010_alphaUp/F");
  myTree->Branch("ZbCHS3010_alphaDown",&nEvent.ZbCHS3010_alphaDown,"ZbCHS3010_alphaDown/F");
  myTree->Branch("ZbCHS1010_alphaDown",&nEvent.ZbCHS1010_alphaDown,"ZbCHS1010_alphaDown/F");
  myTree->Branch("mpf",&nEvent.mpf,"mpf/F");
  myTree->Branch("fake_mpf",&nEvent.fake_mpf,"fake_mpf/F");
  
  myTree->Branch("HasSoftLepton",&nEvent.HasSoftLepton,"HasSoftLepton/O");
  myTree->Branch("SoftLeptonPt",&nEvent.SoftLeptonPt,"SoftLeptonPt/F");
  
  myTree->Branch("dRmerged",&nEvent.dRmerged,"dRmerged/F");
  myTree->Branch("merged",&nEvent.merged,"merged/I");
  myTree->Branch("mergedPt",&nEvent.mergedPt,"mergedPt/F");
  myTree->Branch("mergedEta",&nEvent.mergedEta,"mergedEta/F");
  myTree->Branch("mergedPhoHCalIso2012ConeDR03",&nEvent.mergedPhoHCalIso2012ConeDR03,"mergedPhoHCalIso2012ConeDR03/F");
  myTree->Branch("mergedPhoSigmaIetaIeta",&nEvent.mergedPhoSigmaIetaIeta,"mergedPhoSigmaIetaIeta/F");
  myTree->Branch("mergedpfchargedhadiso",&nEvent.mergedpfchargedhadiso,"mergedpfchargedhadiso/F");
  myTree->Branch("mergedpfneutralhadiso",&nEvent.mergedpfneutralhadiso,"mergedpfneutralhadiso/F");
  myTree->Branch("mergedpfphotoniso",&nEvent.mergedpfphotoniso,"mergedpfphotoniso/F");


  

  counters[EV].setName("Events");
  counters[TR].setName("Triggers");
  counters[MU].setName("Muons");
  counters[EL].setName("Electrons");
  counters[PJ].setName("PFJets");
  counters[PH].setName("Photons");

  // Define counters (so we have them in the right order)
  counters[EV].fill("All events",0.);
  if ( fDataType_ != "mc" ) {
    counters[EV].fill("... pass electron triggers",0.);
    counters[EV].fill("... pass muon triggers",0.);
    counters[EV].fill("... pass EM triggers",0.);
    counters[EV].fill("... pass MET triggers",0.);
  }
  std::string types[3] = { "ee","mm","em" };
  for ( size_t itype=0; itype<3; ++itype ) {
    counters[EV].fill("... "+types[itype]+" pairs",0.); 
    counters[EV].fill("... "+types[itype]+" + 2 jets",0.);
    counters[EV].fill("... "+types[itype]+" + 2 jets + require Z",0.);
    counters[EV].fill("... "+types[itype]+" + 2 jets + require Z + JZB>50",0.);
  }


}


//------------------------------------------------------------------------------
bool momentumComparator(lepton i, lepton j) { return (i.p.Pt()>j.p.Pt()); }

const float JZBAnalysis::GetL5Correction(const int jindex) {
  if(abs(fTR->PFCHSJFlavour[jindex])<4) {
    //dealing with light quark
    return L5corr_qJ->Eval(fTR->PFCHSJPt[jindex]);
  } 
  if(abs(fTR->PFCHSJFlavour[jindex])==4) {
    //dealing c
    return L5corr_cJ->Eval(fTR->PFCHSJPt[jindex]);
  }
  if(abs(fTR->PFCHSJFlavour[jindex])==5) {
    //dealing with b
    return L5corr_bJ->Eval(fTR->PFCHSJPt[jindex]);
  }
  if(abs(fTR->PFCHSJFlavour[jindex])==21) {
    return L5corr_gJ->Eval(fTR->PFCHSJPt[jindex]);
  }
  cerr << "Not able to find any matching L5 function for flavor " << fTR->PFCHSJPt[jindex] << " in event (RLE) " << nEvent.runNum << ":" << nEvent.lumi << ":" << nEvent.eventNum << " (returning correction 1.0) " << endl;
  return 1.0;
}

float JZBAnalysis::GetBWeight(const string WP,const int JetFlavor, const float JetPt, const float JetEta, float &Uncert) {
  Uncert=0.0;
  TH2F *CSV_EfficiencyCorrection;
  TH2F *CSV_EfficiencyCorrectionUncert;
  TH2F *CSV_MisTagCorrection;
  TH2F *CSV_MisTagCorrectionUncert;
  
  bool assigned=false;
  if(WP=="Tight") {
    assigned=true;
    CSV_EfficiencyCorrection       = CSVT_EfficiencyCorrection;
    CSV_EfficiencyCorrectionUncert = CSVT_EfficiencyCorrectionUncert;
    CSV_MisTagCorrection           = CSVT_MisTagCorrection;
    CSV_MisTagCorrectionUncert     = CSVT_MisTagCorrectionUncert;
  }
  if(WP=="Medium") {
    assigned=true;
    CSV_EfficiencyCorrection       = CSVM_EfficiencyCorrection;
    CSV_EfficiencyCorrectionUncert = CSVM_EfficiencyCorrectionUncert;
    CSV_MisTagCorrection           = CSVM_MisTagCorrection;
    CSV_MisTagCorrectionUncert     = CSVM_MisTagCorrectionUncert;
  }
  if(WP=="Loose") {
    assigned=true;
    CSV_EfficiencyCorrection       = CSVL_EfficiencyCorrection;
    CSV_EfficiencyCorrectionUncert = CSVL_EfficiencyCorrectionUncert;
    CSV_MisTagCorrection           = CSVL_MisTagCorrection;
    CSV_MisTagCorrectionUncert     = CSVL_MisTagCorrectionUncert;
  }
  if(!assigned) {
    cout << "You requested the " << WP << " working point. Unfortunately that is not one of the options ... bye bye." << endl;
    assert(assigned);
  }
  
  
    
  
  if(abs(JetFlavor)==4||abs(JetFlavor)==5) {
    //return weight from efficiency correction
    int Bin = CSV_EfficiencyCorrection->FindBin(JetPt,JetEta);
    float weight = CSV_EfficiencyCorrection->GetBinContent(Bin);
    Bin=CSV_EfficiencyCorrectionUncert->FindBin(JetPt,JetEta);
    Uncert=CSV_EfficiencyCorrectionUncert->GetBinContent(Bin);
    if(Uncert<0) Uncert=0.0;
    if(weight<0) {
      Uncert=0.0;
      return 1.0;
    }
    else return weight;
  } else {
    //return weight from mistag correction
    int Bin = CSV_MisTagCorrection->FindBin(JetPt,JetEta);
    float weight = CSV_MisTagCorrection->GetBinContent(Bin);
    Bin = CSV_MisTagCorrectionUncert->FindBin(JetPt,JetEta);
    Uncert=CSV_MisTagCorrectionUncert->GetBinContent(Bin);
    if(Uncert<0) Uncert=0.0;
    if(weight<=0) {
      Uncert=0.0;
      return 1.0;
    } else {
      return weight;
    }
  }
}


//------------------------------------------------------------------------------
vector<lepton> JZBAnalysis::sortLeptonsByPt(vector<lepton>& leptons) {
  
  vector<lepton> theLep = leptons;
  sort (theLep.begin(), theLep.end(), momentumComparator);
  return theLep;  
  
}



//------------------------------------------------------------------------------
//for triggers, check out
// http://fwyzard.web.cern.ch/fwyzard/hlt/summary
const bool JZBAnalysis::passTriggers( std::vector<std::string>& triggerPaths ) {

  bool foundUnprescaled(false);
  bool passed(false);
  for ( size_t i=0; i<triggerPaths.size(); ++i ) {
    if ( GetHLTResult(triggerPaths[i]) ) passed = true;
    if ( GetHLTPrescale(triggerPaths[i]) == 1 ) foundUnprescaled = true;
  }

  // Check if found unprescaled trigger...
//  assert(foundUnprescaled);

  return passed;

}

bool is_neutrino(const int code) {
  if(abs(code)==12) return true; // electron neutrino
  if(abs(code)==14) return true; // muon neutrino
  if(abs(code)==16) return true; // tau neutrino
  if(abs(code)==18) return true; // tau' neutrino
  return false;
}

bool is_charged_lepton(const int code) {
  if(abs(code)==11) return true; // electron
  if(abs(code)==13) return true; // muon
  if(abs(code)==15) return true; // tau
  if(abs(code)==17) return true; // tau'
  return false;
}

//______________________________________________________________________________
void JZBAnalysis::Analyze() {
  // #--- analysis global parameters
  double DRmax=0.4; // veto jets in a cone of DRmax close to the lepton
  counters[EV].fill("All events");
  nEvent.reset();
  // Fill generic information
  nEvent.eventNum  = fTR->Event;
  nEvent.runNum    = fTR->Run;
  nEvent.lumi      = fTR->LumiSection;
  nEvent.totEvents = fTR->GetEntries();

  if(fDataType_ == "mc") // only do this for MC; for data nEvent.reset() has already set both weights to 1 
    {
      if(fisModelScan) {
        nEvent.process=fTR->process;
        nEvent.mGlu=fTR->MassGlu;
        nEvent.mChi=fTR->MassChi;
        nEvent.mLSP=fTR->MassLSP;
        nEvent.mGMSBGlu=fTR->MassChi; // explanation: order in NTuple is wrong for GMSB
        nEvent.mGMSBChi=fTR->MassLSP; // explanation: order in NTuple is wrong for GMSB
        nEvent.mGMSBLSP=fTR->MassGlu; // explanation: order in NTuple is wrong for GMSB
        nEvent.xSMS=fTR->xSMS;
        nEvent.xbarSMS=fTR->xbarSMS;
        nEvent.NPdfs=fTR->NPdfs;
        for(int i=0;i<fTR->NPdfs;i++) nEvent.pdfW[i]=fTR->pdfW[i];
	nEvent.pdfWsum=fTR->pdfWsum; 
      } else {
	//don't attempt to do PURW for model scans
	nEvent.PUweight     = GetPUWeight(fTR->PUnumTrueInteractions);
	nEvent.PUweightUP   = GetPUWeightUp(fTR->PUnumTrueInteractions);
	nEvent.PUweightDOWN = GetPUWeightDown(fTR->PUnumTrueInteractions);
	nEvent.weight     = nEvent.PUweight;
	weight_histo->Fill(1,nEvent.PUweight);
      }
      
      nEvent.EventFlavor=DetermineFlavor(fdoGenInfo,fTR);
      nEvent.EventZToTaus=DecaysToTaus(fdoGenInfo,fTR);
      
     // the following part makes sense for all MC - not only for scans (though for scans imposedx/realx make more sense)
	float chimass=0;
	int nchimass=0;
	float lspmass=0;
	int nlspmass=0;
	float glumass=0;
	int nglumass=0;
	
	float genZpt=0,genZeta=0,genZphi=0,genZM=0;
	float genZ2pt=0,genZ2eta=0,genZ2phi=0,genZ2M=0;
	int Zprompt1=0,Zprompt2=0;
	int Promptness[5];
	vector<TLorentzVector> LSPvecs;
	vector<TLorentzVector> LSPMothervecs;
	vector<int> LSPMother;
	vector<float> LSPMotherPt;

	TLorentzVector summedLSPs;
	int nGenParticles=fTR->nGenParticles;
	if(nGenParticles<2||nGenParticles>2000) {
		//this happens if you use an old file or one that doesn't contain the necessary generator information.
		if(fdoGenInfo) cerr << "WATCH OUT : GENERATOR INFORMATION HAS BEEN DISABLED BECAUSE THE NUMBER OF GEN PARTICLES WAS TOO LOW (" << nGenParticles << ")" << endl;
		fdoGenInfo=false;
		nGenParticles=0;
	}


	ContainsSoftLepton();


	bool wecare=false;
	for(int i=0;i<nGenParticles&&fdoGenInfo;i++) {
	  int thisParticleId = fTR->genInfoId[i];
	  if(fTR->genInfoStatus[i]!=3) continue;
	  if(fdoGenInfo&&abs(thisParticleId)==23) {
	    //dealing with a Z
	    int motherIndex=fTR->genInfoMo1[i];
	    if(motherIndex>=0) nEvent.SourceOfZ[nEvent.nZ]=fTR->genInfoId[motherIndex];
	    nEvent.nZ++;
	    for(int da=i+1;da<fTR->nGenParticles;da++) {
	      if(fTR->genInfoMo1[da]==i) {
		//dealing with a daughter
		if(abs(fTR->genInfoId[da])<10) nEvent.DecayCode+=100;
		if(is_neutrino(abs(fTR->genInfoId[da]))) nEvent.DecayCode+=10;
		if(is_charged_lepton(abs(fTR->genInfoId[da]))) {
		  nEvent.DecayCode+=1;
		  if(fTR->genInfoPt[i]>genZpt) {
		    genZ2pt=genZpt;
		    genZ2M=genZM;
		    genZ2eta=genZeta;
		    genZ2phi=genZphi;
		    Zprompt2=Zprompt1;
		    
		    genZpt=fTR->genInfoPt[i];
		    genZM=fTR->genInfoM[i];
		    genZeta=fTR->genInfoEta[i];
		    genZphi=fTR->genInfoPhi[i];
		    Zprompt1=fTR->PromptnessLevel[i];
		  } else {
		    if(fTR->genInfoPt[i]>genZ2pt) {
		      genZ2pt=fTR->genInfoPt[i];
		      genZ2M=fTR->genInfoM[i];
		      genZ2eta=fTR->genInfoEta[i];
		      genZ2phi=fTR->genInfoPhi[i];
		      Zprompt2=fTR->PromptnessLevel[i];
		    }
		  }//end of if fTR->genInfoPt[i]>genZpt)
		}//end of if leptonic decay
	      }//end of if daughter
	    }//end of daughter search
	  }//end of Z case
	  
	  if(abs(thisParticleId)==1000021) {//mglu
	    glumass+=fTR->genInfoM[i];
	    nglumass++;
	  }
	  if(abs(thisParticleId)==1000022) {//mlsp
	    LSPMother.push_back(fTR->genInfoMo1[i]);
	    LSPMotherPt.push_back(fTR->genInfoPt[fTR->genInfoMo1[i]]);
	    lspmass+=fTR->genInfoM[i];
	    nlspmass++;
	    TLorentzVector newLSP;
	    newLSP.SetPtEtaPhiM(fTR->genInfoPt[i],fTR->genInfoEta[i],fTR->genInfoPhi[i],fTR->genInfoM[i]);
	    LSPvecs.push_back(newLSP);
	    if(LSPvecs.size()==1) summedLSPs.SetPtEtaPhiM(fTR->genInfoPt[i],fTR->genInfoEta[i],fTR->genInfoPhi[i],fTR->genInfoM[i]);
	    else summedLSPs=summedLSPs+newLSP;
	    Promptness[nEvent.nLSPs]=fTR->PromptnessLevel[i];
	    nEvent.nLSPs++;
	  }
	  if(abs(thisParticleId)==1000023) {//mchi
	    chimass+=fTR->genInfoM[i];
	    nchimass++;
	    TLorentzVector thismom;
	    thismom.SetPtEtaPhiM(fTR->genInfoPt[i],fTR->genInfoEta[i],fTR->genInfoPhi[i],fTR->genInfoM[i]);
	    LSPMothervecs.push_back(thismom);
	  }
	}// done with gen info loop

	if(fdoGenInfo&&fisModelScan) {
		TLorentzVector pureGenMETvector(fTR->GenMETpx,fTR->GenMETpy,0,0);

		if(nEvent.nLSPs==2) nEvent.angleLSPLSP=LSPvecs[0].Angle(LSPvecs[1].Vect());
		if(nEvent.nLSPs==2) nEvent.angleLSPLSP2d=LSPvecs[0].DeltaPhi(LSPvecs[1]);
		
		nEvent.dphiSumLSPgenMET=summedLSPs.DeltaPhi(pureGenMETvector);
		nEvent.SumLSPEta=summedLSPs.Eta();
		nEvent.SumLSPPhi=summedLSPs.Phi();
		nEvent.absvalSumLSP=summedLSPs.Pt();

		TLorentzVector pureGenZvector;
		pureGenZvector.SetPtEtaPhiM(genZpt,genZeta,genZphi,genZM);

		
		float LSPdRa = LSPvecs[0].DeltaR(pureGenZvector);
		float LSPdRb = LSPvecs[1].DeltaR(pureGenZvector);
		nEvent.ZPromptnessLevel[0]=Zprompt1;
		nEvent.ZPromptnessLevel[1]=Zprompt2;
		if(LSPdRa<LSPdRb) {
			nEvent.LSPPromptnessLevel[0]=Promptness[0];
			nEvent.LSPPromptnessLevel[1]=Promptness[1];
			nEvent.LSP1pt=LSPvecs[0].Pt();
			nEvent.LSP2pt=LSPvecs[1].Pt();
			nEvent.LSP1Mo=LSPMother[0];
			nEvent.LSP2Mo=LSPMother[1];
			nEvent.LSP1Mopt=LSPMotherPt[0];
			nEvent.LSP2Mopt=LSPMotherPt[1];
			nEvent.angleLSPZ=LSPvecs[0].Angle(pureGenZvector.Vect());
			nEvent.angleLSPZ2d=LSPvecs[0].DeltaPhi(pureGenZvector);
			if(abs(LSPMotherPt[0]-LSPMothervecs[0].Pt())<abs(LSPMotherPt[0]-LSPMothervecs[1].Pt())) {
			  nEvent.angleChi2Z2d=LSPMothervecs[0].DeltaPhi(LSPvecs[0]);
			  nEvent.angleChi2Z=LSPMothervecs[0].Angle(LSPvecs[0].Vect());
			} else {
			  nEvent.angleChi2Z2d=LSPMothervecs[1].DeltaPhi(LSPvecs[0]);
			  nEvent.angleChi2Z=LSPMothervecs[1].Angle(LSPvecs[0].Vect());
			}
		} else {
			nEvent.LSPPromptnessLevel[0]=Promptness[1];
			nEvent.LSPPromptnessLevel[1]=Promptness[0];
			nEvent.LSP1pt=LSPvecs[1].Pt();
			nEvent.LSP2pt=LSPvecs[0].Pt();
			nEvent.LSP1Mo=LSPMother[1];
			nEvent.LSP2Mo=LSPMother[0];
			nEvent.LSP1Mopt=LSPMotherPt[1];
			nEvent.LSP2Mopt=LSPMotherPt[0];
			nEvent.angleLSPZ=LSPvecs[1].Angle(pureGenZvector.Vect());
			nEvent.angleLSPZ2d=LSPvecs[1].DeltaPhi(pureGenZvector);
			if(abs(LSPMotherPt[1]-LSPMothervecs[0].Pt())<abs(LSPMotherPt[1]-LSPMothervecs[1].Pt())) {
			  nEvent.angleChi2Z2d=LSPMothervecs[0].DeltaPhi(LSPvecs[1]);
			  nEvent.angleChi2Z=LSPMothervecs[0].Angle(LSPvecs[1].Vect());
			} else {
			  nEvent.angleChi2Z2d=LSPMothervecs[1].DeltaPhi(LSPvecs[1]);
			  nEvent.angleChi2Z=LSPMothervecs[1].Angle(LSPvecs[1].Vect());
			}
		}

		TLorentzVector pureGenZ2vector;
		pureGenZ2vector.SetPtEtaPhiM(genZ2pt,genZ2eta,genZ2phi,genZ2M);
		nEvent.pureGeneratorJZB=(-pureGenMETvector-pureGenZvector).Pt() - pureGenZvector.Pt();
		nEvent.pureGeneratorMet=pureGenMETvector.Pt();
		nEvent.pureGeneratorMetPhi=pureGenMETvector.Phi();
		nEvent.pureGeneratorSumJetPt=(-pureGenMETvector-pureGenZvector).Pt();
		nEvent.pureGeneratorSumJetEta=(-pureGenMETvector-pureGenZvector).Eta();
		nEvent.pureGeneratorSumJetPhi=(-pureGenMETvector-pureGenZvector).Phi();


		nEvent.pureGeneratorZpt=genZpt;
		nEvent.pureGeneratorZM=genZM;
		nEvent.pureGeneratorZphi=genZphi;
		nEvent.pureGeneratorZeta=genZeta;
		nEvent.pure2ndGeneratorJZB=(-pureGenMETvector-pureGenZ2vector).Pt() - pureGenZ2vector.Pt();
		nEvent.pure2ndGeneratorZpt=pureGenZ2vector.Pt();

		if(genZpt<0.01) nEvent.pureGeneratorJZB=0; // in case there is no leptonic Z
	}//end of if(fdoGenInfo&&fisModelScan)

	
	if(nchimass>0&&nlspmass>0&&nglumass>0)  nEvent.realx=(chimass/nchimass - lspmass/nlspmass)/(glumass/nglumass-lspmass/nlspmass);
	//at this point we use the fact that one of the three bits of information in the LHE event comment is the imposed x - the only question is which one ;-) 
	//note: the bit of information in the comment is actually xbar, so we need to store 1-xbar to get our definition of x. 
	if(nEvent.mGlu>0 && nEvent.mGlu<1) nEvent.imposedx=1-nEvent.mGlu;
	if(nEvent.mChi>0 && nEvent.mChi<1) nEvent.imposedx=1-nEvent.mChi;
	if(nEvent.mLSP>0 && nEvent.mLSP<1) nEvent.imposedx=1-nEvent.mLSP;
  } // end of mc if
  
  // Trigger information
  nEvent.passed_triggers=0;
  if ( fDataType_ != "mc" ) nEvent.is_data=true;
  if ( passTriggers(elTriggerPaths) ) 
    {
      counters[EV].fill("... pass electron triggers");
      nEvent.passed_triggers=1;
      nEvent.trigger_bit |= 1;
    } 
  if ( passTriggers(muTriggerPaths) )
    {
      counters[EV].fill("... pass muon triggers");
      nEvent.passed_triggers=1;
      nEvent.trigger_bit |= (1<<1);
    } 
  if ( passTriggers(emTriggerPaths) )
    {
      counters[EV].fill("... pass EM triggers");
      nEvent.passed_triggers=1;
      nEvent.trigger_bit |= (1<<2);
    }
  if ( passTriggers(meTriggerPaths) )
    {
      counters[EV].fill("... pass ME triggers");
      nEvent.passed_triggers=1;
      nEvent.trigger_bit |= (1<<3);
    }
  if ( passTriggers(metTriggerPaths) )
    {
      counters[EV].fill("... pass MET triggers");
      nEvent.passed_triggers=1;
      nEvent.trigger_bit |= (1<<4);
    }
  if ( passTriggers(htTriggerPaths) )
    {
      counters[EV].fill("... pass HT triggers");
      nEvent.passed_triggers=1;
      nEvent.trigger_bit |= (1<<5);
    }
    
  // Event filter information
  nEvent.passed_filters = passFilters( nEvent.filter_bit );

  // Check if we find an OSSF pair in the acceptance (and if it is coming from a Z)
  bool isMC = (fDataType_ == "mc");
  if ( isMC ) { GeneratorInfo(); }
    
  // #--- Vertex info
  nEvent.numVtx = fTR->NVrtx;
  float rho = sqrt(fTR->PrimVtxx*fTR->PrimVtxx + fTR->PrimVtxy*fTR->PrimVtxy);
  if(fTR->PrimVtxGood) nEvent.goodVtx |=2; // save bits of vertex quality
  if (   fTR->PrimVtxGood==0 && fTR->PrimVtxIsFake==0 
         && fTR->PrimVtxNdof>4  && fabs(fTR->PrimVtxz)<24 && rho<2)
    nEvent.goodVtx |=4;
  
  // Good event requirement: essentially vertex requirements
  if ( !IsGoodEvent() ) {
    if (isMC&&!fmakeSmall) myTree->Fill();
    return;
  }
  counters[EV].fill("... pass good event requirements");

  vector<lepton> leptons;

  TLorentzVector genZvector; // To store the true Z vector
  vector<lepton> gLeptons;   // lepton collection corrected for gammas

  // #--- muon loop
  for(int muIndex=0;muIndex<fTR->NMus;muIndex++)
    {
      counters[MU].fill("All mus");
      if(IsCustomMu2012(muIndex))
        {
          counters[MU].fill("... pass mu selection");
          float px= fTR->MuPx[muIndex];
          float py= fTR->MuPy[muIndex];
          float pz= fTR->MuPz[muIndex];
          float energy =  fTR->MuE[muIndex];
          TLorentzVector tmpVector(px,py,pz,energy);
          int tmpCharge = fTR->MuCharge[muIndex];
          float muonIso = MuPFIso(muIndex);
          lepton tmpLepton;
          tmpLepton.p = tmpVector;
          tmpLepton.charge = tmpCharge;
          tmpLepton.index = muIndex;
          tmpLepton.iso   = muonIso;
          tmpLepton.type = 1;
          tmpLepton.genPt = 0.;
          tmpLepton.ElCInfoIsGsfCtfCons=true;
          tmpLepton.ElCInfoIsGsfCtfScPixCons=true;
          tmpLepton.ElCInfoIsGsfScPixCons=true;
          leptons.push_back(tmpLepton);
        }
    }

  // #--- electron loop
  for(int elIndex=0;elIndex<fTR->NEles;elIndex++)
    {
      counters[EL].fill("All eles");
      if(IsCustomEl2012(elIndex))	
        {
          counters[EL].fill("... pass e selection");
          float px= fTR->ElPx[elIndex];
          float py= fTR->ElPy[elIndex];
          float pz= fTR->ElPz[elIndex];
          float energy =  fTR->ElE[elIndex];
          TLorentzVector tmpVector(px,py,pz,energy);
          int tmpCharge=fTR->ElCharge[elIndex];
          double pedestal=0.;
          if ( fabs(fTR->ElEta[elIndex]) < 1.479 ) pedestal = 1.0;
          double pfIso = ElPFIso(elIndex);
          lepton tmpLepton;
          tmpLepton.p = tmpVector;
          tmpLepton.charge = tmpCharge;
          tmpLepton.index = elIndex;
          tmpLepton.iso = pfIso;
          tmpLepton.type = 0;
          tmpLepton.genPt = 0.;
          tmpLepton.ElCInfoIsGsfCtfCons=fTR->ElCInfoIsGsfCtfCons[elIndex];
          tmpLepton.ElCInfoIsGsfCtfScPixCons=fTR->ElCInfoIsGsfCtfScPixCons[elIndex];
          tmpLepton.ElCInfoIsGsfScPixCons=fTR->ElCInfoIsGsfScPixCons[elIndex];
          leptons.push_back(tmpLepton);
        }
    }

  // Sort the leptons by Pt and select the two opposite-signed ones with highest Pt
  vector<lepton> sortedGoodLeptons = sortLeptonsByPt(leptons);

  if(sortedGoodLeptons.size() < 2) {
    if (isMC&&!fmakeSmall) myTree->Fill();
    return;
  }

    
  counters[EV].fill("... has at least 2 leptons");
  int PosLepton1 = 0;
  int PosLepton2 = 1;
    
  int TriLepton1 = 0;
  int TriLepton2 = 0;
  int TriLepton3 = 0;
  
  int BadTriLepton1 = 0;
  int BadTriLepton2 = 0;
  int BadTriLepton3 = 0;
  
//  // Check for SS combination
//  for(; PosLepton2 < sortedGoodgLeptons.size(); PosLepton2++) {
//    if(sortedGoodgLeptons[0].charge*sortedGoodgLeptons[PosLepton2].charge>0) break;
//  }
  
  // Check for OS combination
  for(; PosLepton2 < sortedGoodLeptons.size(); PosLepton2++) {
    if(sortedGoodLeptons[0].charge*sortedGoodLeptons[PosLepton2].charge<0) break;
  }
  if(PosLepton2 == sortedGoodLeptons.size()) {
    if (isMC&&!fmakeSmall) myTree->Fill();
    return;
  }
  
  counters[EV].fill("... has at least 2 OS leptons");
  
  float mindiff=-1;
  int BestCandidate=sortedGoodLeptons.size();
  
  // Trileptons: the highest pt lepton that fired the trigger is our first lepton!
  nEvent.tri_MatchFound=false;
  
  for(int ilep=0;ilep<sortedGoodLeptons.size();ilep++) {
    if(MatchTrigger(&sortedGoodLeptons[ilep])) {
      TriLepton1=ilep;
      nEvent.tri_MatchFound=true;
      break;
    }
  }
  
  if(!nEvent.tri_MatchFound) TriLepton1=0;
  BadTriLepton1=TriLepton1;
  
  for(; TriLepton2 < sortedGoodLeptons.size(); TriLepton2++) {
    if(TriLepton2==TriLepton1) continue;
    if(sortedGoodLeptons[0].charge*sortedGoodLeptons[TriLepton2].charge<0 &&  sortedGoodLeptons[0].type==sortedGoodLeptons[TriLepton2].type ) {
      float curr_mindiff=abs((sortedGoodLeptons[0].p+sortedGoodLeptons[TriLepton2].p).M()-91.2);
      if(mindiff<0||curr_mindiff<mindiff) {
	mindiff=curr_mindiff;
	BestCandidate=TriLepton2;
      }
    }
  }
  TriLepton2=BestCandidate;
  
  for(; TriLepton3 < sortedGoodLeptons.size(); TriLepton3++) {
    if(TriLepton3==TriLepton1 || TriLepton3==TriLepton2) continue;
    else break; // pick the first lepton that's left!
  }
  
  mindiff=-1;
  BestCandidate=sortedGoodLeptons.size();
  
  // Trileptons
  for(; BadTriLepton2 < sortedGoodLeptons.size(); BadTriLepton2++) {
    if(BadTriLepton2==BadTriLepton1) continue;
    if(sortedGoodLeptons[0].charge*sortedGoodLeptons[BadTriLepton2].charge<0 ) { //ignore flavor to get mismatch rate
      float curr_mindiff=abs((sortedGoodLeptons[0].p+sortedGoodLeptons[BadTriLepton2].p).M()-91.2);
      if(mindiff<0||curr_mindiff<mindiff) {
	mindiff=curr_mindiff;
	BestCandidate=BadTriLepton2;
      }
    }
  }
  BadTriLepton2=BestCandidate;
  
  for(; BadTriLepton3 < sortedGoodLeptons.size(); BadTriLepton3++) {
    if(BadTriLepton3==BadTriLepton1 || BadTriLepton3==BadTriLepton2) continue;
    else break; // pick the first lepton that's left!
  }
  

  
  
  //***************//
  
  // Preselection
  if(sortedGoodLeptons[PosLepton1].p.Pt() >= firstLeptonPtCut && sortedGoodLeptons[PosLepton2].p.Pt() >= secondLeptonPtCut) {
    nEvent.eta1 = sortedGoodLeptons[PosLepton1].p.Eta();
    nEvent.pt1 = sortedGoodLeptons[PosLepton1].p.Pt();
    nEvent.iso1 = sortedGoodLeptons[PosLepton1].iso;
    nEvent.phi1 = sortedGoodLeptons[PosLepton1].p.Phi();
    nEvent.ch1 = sortedGoodLeptons[PosLepton1].charge;
    nEvent.id1 = sortedGoodLeptons[PosLepton1].type; //??????
    nEvent.chid1 = (sortedGoodLeptons[PosLepton1].type+1)*sortedGoodLeptons[PosLepton1].charge;
//    nEvent.isConv1 = IsConvertedPhoton(sortedGoodLeptons[PosLepton1].index);
      
    nEvent.eta2 = sortedGoodLeptons[PosLepton2].p.Eta();
    nEvent.pt2 = sortedGoodLeptons[PosLepton2].p.Pt();
    nEvent.iso2 = sortedGoodLeptons[PosLepton2].iso;
    nEvent.phi2 = sortedGoodLeptons[PosLepton2].p.Phi();
    nEvent.ch2 = sortedGoodLeptons[PosLepton2].charge;
    nEvent.id2 = sortedGoodLeptons[PosLepton2].type; //??????
    nEvent.chid2 = (sortedGoodLeptons[PosLepton2].type+1)*sortedGoodLeptons[PosLepton2].charge;
//    nEvent.isConv2 = IsConvertedPhoton(sortedGoodLeptons[PosLepton2].index);
    
    nEvent.l1l2dR=sortedGoodLeptons[PosLepton1].p.DeltaR(sortedGoodLeptons[PosLepton2].p);
    
    nEvent.minc=sortedGoodLeptons[PosLepton2].p.Pt()+sortedGoodLeptons[PosLepton1].p.Pt();
    nEvent.mll=(sortedGoodLeptons[PosLepton2].p+sortedGoodLeptons[PosLepton1].p).M();
    nEvent.phi=(sortedGoodLeptons[PosLepton2].p+sortedGoodLeptons[PosLepton1].p).Phi();
    nEvent.eta=(sortedGoodLeptons[PosLepton2].p+sortedGoodLeptons[PosLepton1].p).Eta();
    nEvent.E=(sortedGoodLeptons[PosLepton2].p+sortedGoodLeptons[PosLepton1].p).E();
    nEvent.pt=(sortedGoodLeptons[PosLepton2].p+sortedGoodLeptons[PosLepton1].p).Pt();
    nEvent.dphi=sortedGoodLeptons[PosLepton2].p.DeltaPhi(sortedGoodLeptons[PosLepton1].p);
    
    nEvent.ElCInfoIsGsfCtfCons=sortedGoodLeptons[PosLepton2].ElCInfoIsGsfCtfCons&&sortedGoodLeptons[PosLepton1].ElCInfoIsGsfCtfCons;
    nEvent.ElCInfoIsGsfCtfScPixCons=sortedGoodLeptons[PosLepton2].ElCInfoIsGsfCtfScPixCons&&sortedGoodLeptons[PosLepton1].ElCInfoIsGsfCtfScPixCons;
    nEvent.ElCInfoIsGsfScPixCons=sortedGoodLeptons[PosLepton2].ElCInfoIsGsfScPixCons&&sortedGoodLeptons[PosLepton1].ElCInfoIsGsfScPixCons;

    float lepweightErr;
    float lepweight=GetLeptonWeight(nEvent.id1,nEvent.pt1,nEvent.eta1,nEvent.id2,nEvent.pt2,nEvent.eta2,lepweightErr);
    
    bool softMuon = false;
    for(int muIndex=0;muIndex<fTR->NMus;muIndex++) {
      if(IsSoftMuon(muIndex)) {
	softMuon = true;
	break;
      }
    }
    nEvent.softMuon = softMuon;
    if (isMC) {
      bool softMuonMC = false;
      for(int muIndex=0;muIndex<fTR->NMus;muIndex++) {
	if(fabs(fTR->MuGenMID[muIndex]) == 24 && fabs(fTR->MuGenGMID[muIndex]) == 5) {
	  softMuonMC = true;
	  break;
	}
      }
      nEvent.softMuonMC = softMuonMC;
    }
    
    if (isMC) {
//      nEvent.weight=nEvent.weight*lepweight;
      nEvent.weightEffDown=nEvent.weight*(lepweight-lepweightErr);
      nEvent.weightEffUp=nEvent.weight*(lepweight+lepweightErr);
      nEvent.Efficiencyweightonly=lepweight;
    }

  } else {
    //If there are less than two leptons the event is not considered
    if (isMC&&!fmakeSmall) myTree->Fill();
    return;
  }
  counters[EV].fill("... pass dilepton pt selection");
        
  // #--- construct different recoil models, initial the recoil vector will hold only the sum over the hard jets, only in the end we will add-up the lepton system

  // --- construct met vectors here
  float pfMETpx = fTR->PFMETpx;
  float pfMETpy = fTR->PFMETpy;
  
  float tcMETpx = fTR->TCMETpx;
  float tcMETpy = fTR->TCMETpy;
  
  float type1METpx = fTR->PFType1METpx;
  float type1METpy = fTR->PFType1METpy;
  
  //Use this factor to stretch MET and, consequently, also JZB.
  nEvent.fact = 1.0;
  if(isMC && fIsDY && ! nEvent.EventZToTaus) nEvent.fact = r->Gaus(1.08,0.1); // only apply this to Z->ee and Z->mm (i.e. only a PART of DY)

  nEvent.minc+= fTR->PFMET;

  TLorentzVector pfMETvector(nEvent.fact*pfMETpx,nEvent.fact*pfMETpy,0,0);
  TLorentzVector tcMETvector(nEvent.fact*tcMETpx,nEvent.fact*tcMETpy,0,0);
  TLorentzVector type1METvector(nEvent.fact*type1METpx,nEvent.fact*type1METpy,0,0);
  TLorentzVector Cleantype1METvector(type1METpx,type1METpy,0,0);
  TLorentzVector sumOfPFJets(0,0,0,0);
  TLorentzVector zVector;
  zVector.SetPtEtaPhiE(nEvent.pt,nEvent.eta,nEvent.phi,nEvent.E);

  nEvent.ZbCHS3010_pfJetGoodNum=0;
  nEvent.ZbCHS1010_pfJetGoodNum=0;
  nEvent.ZbCHS3010_pfJetGoodNumBtag=0;
    // #--- PF jet loop (with CHS - not what we usually use)
  vector<lepton> pfCHSGoodJets;
  for(int i =0 ; i<fTR->PFCHSNJets;i++) // PF jet loop
    {
      if(i==jMax){cout<<"max Num of jets (CHS) was reached"<<endl; break;}
      float jpt = fTR->PFCHSJPt[i];
      float jeta = fTR->PFCHSJEta[i];
      float jphi = fTR->PFCHSJPhi[i];
      float jenergy = fTR->PFCHSJE[i];
      bool isJetID = fTR->PFCHSJIDLoose[i];
      
      float smeared_jpt = smearedJetPt(jpt,jeta,jphi);
      TLorentzVector aJet(0,0,0,0);
      aJet.SetPtEtaPhiE(jpt, jeta, jphi, jenergy);

      // lepton-jet cleaning
      if ( fFullCleaning_ ) {
        // Remove jet close to any lepton
        bool isClean(true);
        for ( size_t ilep = 0; ilep<sortedGoodLeptons.size(); ++ilep )
          if ( aJet.DeltaR(sortedGoodLeptons[ilep].p)<DRmax) isClean=false;
        if ( !isClean ) continue;
      } else {
        // Remove jet close to leptons from Z candidate
        if(aJet.DeltaR(sortedGoodLeptons[PosLepton1].p)<DRmax) continue;
        if(aJet.DeltaR(sortedGoodLeptons[PosLepton2].p)<DRmax) continue;
      }

      if ( !(fabs(jeta)<3.0 ) ) continue;

      if( (nEvent.ZbCHS3010_pfJetGoodNum==0 && jpt>30 && isJetID && abs(jeta)<2.4) || (nEvent.ZbCHS3010_pfJetGoodNum>0 && jpt>10 && isJetID && abs(jeta)<2.4)) {
	//Z+b selection with 30 GeV leading jet, 10 GeV sub-leading jet
	if(nEvent.ZbCHS3010_pfJetGoodNum==0 && isMC) {
	  float Uncert;
	  nEvent.ZbCHS3010_BTagWgtT     = GetBWeight("Tight",abs(fTR->PFCHSJFlavour[i]), jpt, abs(jeta),Uncert);
	  nEvent.ZbCHS3010_BTagWgtTDown = nEvent.ZbCHS3010_BTagWgtT-Uncert;
	  nEvent.ZbCHS3010_BTagWgtTUp   = nEvent.ZbCHS3010_BTagWgtT+Uncert;
	  nEvent.ZbCHS3010_BTagWgtM     = GetBWeight("Medium",abs(fTR->PFCHSJFlavour[i]), jpt, abs(jeta),Uncert);
	  nEvent.ZbCHS3010_BTagWgtMDown = nEvent.ZbCHS3010_BTagWgtM-Uncert;
	  nEvent.ZbCHS3010_BTagWgtMUp   = nEvent.ZbCHS3010_BTagWgtM+Uncert;
	  nEvent.ZbCHS3010_BTagWgtL     = GetBWeight("Loose",abs(fTR->PFCHSJFlavour[i]), jpt, abs(jeta),Uncert);
	  nEvent.ZbCHS3010_BTagWgtLDown = nEvent.ZbCHS3010_BTagWgtL-Uncert;
	  nEvent.ZbCHS3010_BTagWgtLUp   = nEvent.ZbCHS3010_BTagWgtL+Uncert;
	  
	  nEvent.ZbCHS3010_LeadingJetIsPu=IsPUJet(jpt,jeta,jphi);
	}

	if(nEvent.ZbCHS3010_pfJetGoodNum==1) {
	  nEvent.ZbCHS3010_SubLeadingJetIsPu=IsPUJet(jpt,jeta,jphi);
	  if(smeared_jpt>0) {
	    nEvent.ZbCHS3010_alphaUp = smeared_jpt;
	    nEvent.ZbCHS3010_alphaDown = jpt*jpt/smeared_jpt;
	  } else {
	    //jet did not correspond to any gen jet - smearing doesn't make much sense.
	    nEvent.ZbCHS3010_alphaUp = 10e7;
	    nEvent.ZbCHS3010_alphaDown = 10e7;
	  }
	}
	nEvent.ZbCHS3010_bTagProbCSVBP[nEvent.ZbCHS3010_pfJetGoodNum]=fTR->JnewPFCombinedSecondaryVertexBPFJetTags[i];
	if(nEvent.ZbCHS3010_bTagProbCSVBP[nEvent.ZbCHS3010_pfJetGoodNum]>0.244) {
	  nEvent.ZbCHS3010_pfBJetDphiZ[nEvent.ZbCHS3010_pfJetGoodNumBtag]=aJet.DeltaPhi(zVector);
	  nEvent.ZbCHS3010_pfJetGoodNumBtag++;
	}
	nEvent.ZbCHS3010_pfJetGoodEta[nEvent.ZbCHS3010_pfJetGoodNum]=jeta;
	nEvent.ZbCHS3010_pfJetDphiZ[nEvent.ZbCHS3010_pfJetGoodNum]=aJet.DeltaPhi(zVector);
	nEvent.ZbCHS3010_pfJetGoodPt[nEvent.ZbCHS3010_pfJetGoodNum]=jpt;
	nEvent.ZbCHS3010_L5corr=1.0;
	if(isMC) nEvent.ZbCHS3010_L5corr = GetL5Correction(i);
	nEvent.ZbCHS3010_pfJetGoodNum++;
      }

      if(jpt>10 && isJetID && abs(jeta)<2.4) {
	//Z+b selection with 10 GeV leading jet, 10 GeV sub-leading jet
	if(nEvent.ZbCHS1010_pfJetGoodNum==0 && isMC) {
	  float Uncert;
	  nEvent.ZbCHS1010_BTagWgtT     = GetBWeight("Tight",abs(fTR->PFCHSJFlavour[i]), jpt, abs(jeta),Uncert);
	  nEvent.ZbCHS1010_BTagWgtTDown = nEvent.ZbCHS1010_BTagWgtT - Uncert;
	  nEvent.ZbCHS1010_BTagWgtTUp   = nEvent.ZbCHS1010_BTagWgtT + Uncert;
	  nEvent.ZbCHS1010_BTagWgtM     = GetBWeight("Medium",abs(fTR->PFCHSJFlavour[i]), jpt, abs(jeta),Uncert);
	  nEvent.ZbCHS1010_BTagWgtMDown = nEvent.ZbCHS1010_BTagWgtM - Uncert;
	  nEvent.ZbCHS1010_BTagWgtMUp   = nEvent.ZbCHS1010_BTagWgtM + Uncert;
	  nEvent.ZbCHS1010_BTagWgtL     = GetBWeight("Loose",abs(fTR->PFCHSJFlavour[i]), jpt, abs(jeta),Uncert);
	  nEvent.ZbCHS1010_BTagWgtLDown = nEvent.ZbCHS1010_BTagWgtL - Uncert;
	  nEvent.ZbCHS1010_BTagWgtLUp   = nEvent.ZbCHS1010_BTagWgtL + Uncert;
	  
	  nEvent.ZbCHS1010_LeadingJetIsPu=IsPUJet(jpt,jeta,jphi);
	}
	  
	if(nEvent.ZbCHS1010_pfJetGoodNum==1) {
	  nEvent.ZbCHS1010_SubLeadingJetIsPu=IsPUJet(jpt,jeta,jphi);
	  if(smeared_jpt>0) {
	    nEvent.ZbCHS1010_alphaUp = smeared_jpt;
	    nEvent.ZbCHS1010_alphaDown = jpt*jpt/smeared_jpt;
	  } else {
	    //jet did not correspond to any gen jet - smearing doesn't make much sense.
	    nEvent.ZbCHS1010_alphaUp = 10e7;
	    nEvent.ZbCHS1010_alphaDown = 10e7;
	  }
	}
	
	nEvent.ZbCHS1010_bTagProbCSVBP[nEvent.ZbCHS1010_pfJetGoodNum]=fTR->JnewPFCombinedSecondaryVertexBPFJetTags[i];
	if(nEvent.ZbCHS1010_bTagProbCSVBP[nEvent.ZbCHS1010_pfJetGoodNum]>0.244) {
	  nEvent.ZbCHS1010_pfBJetDphiZ[nEvent.ZbCHS1010_pfJetGoodNumBtag]=aJet.DeltaPhi(zVector);
	  nEvent.ZbCHS1010_pfJetGoodNumBtag++;
	}
	nEvent.ZbCHS1010_pfJetGoodEta[nEvent.ZbCHS1010_pfJetGoodNum]=jeta;
	nEvent.ZbCHS1010_pfJetDphiZ[nEvent.ZbCHS1010_pfJetGoodNum]=aJet.DeltaPhi(zVector);
	nEvent.ZbCHS1010_pfJetGoodPt[nEvent.ZbCHS1010_pfJetGoodNum]=jpt;
	nEvent.ZbCHS1010_pfJetGoodNum++;
      }
      
      //JZB case!!!
      if(jpt>40 && isJetID) nEvent.pfJetGoodNum40CHS++;
      
    }

  nEvent.pfJetNum=0;
  nEvent.pfJetGoodNum30=0;
  nEvent.pfJetGoodNumBtag30=0;
  nEvent.pfJetGoodNumBtag40=0;
  nEvent.pfJetGoodNumBtag30_Tight=0;
  nEvent.pfJetGoodNumBtag40_Tight=0;
  nEvent.pfJetGoodNumBtag30_Loose=0;
  nEvent.pfJetGoodNumBtag40_Loose=0;
  nEvent.pfJetGoodNum40=0;
  nEvent.pfJetGoodNum50=0;
  nEvent.pfJetGoodNum60=0;
  
  TLorentzVector ResidualMet(pfMETpx,pfMETpy,0,0);
  TLorentzVector VariedMet(0,0,0,0);
  TLorentzVector Type1ResidualMet(type1METpx,type1METpy,0,0);
  TLorentzVector Type1VariedMet(0,0,0,0);
  
  // #--- PF jet loop (this is what we use)
  vector<lepton> pfGoodJets;
  for(int i =0 ; i<fTR->NJets;i++) // PF jet loop
    {
      counters[PJ].fill("All PF jets");
      if(i==jMax){cout<<"max Num was reached"<<endl; break;}
	
      float jpt = fTR->JPt[i];
      float jeta = fTR->JEta[i];
      float jphi = fTR->JPhi[i];
      float jpx = fTR->JPx[i];
      float jpy = fTR->JPy[i];
      float jpz = fTR->JPz[i];
      float jenergy = fTR->JE[i];
      float jesC = fTR->JEcorr[i];
      bool  isJetID = IsGoodBasicPFJet(i,0.0,5.0);
      int   ntracksch = fTR->JNAssoTracks[i];  
      int   ntracksne = fTR->JNConstituents[i];


//      jpt=jpt/jesC;
//      jenergy=jenergy/jesC;
//      jpx/=jesC;
//      jpy/=jesC;
//      jpz/=jesC;
//      fJetCorrector->setJetEta(jeta);
//      fJetCorrector->setJetPt(jpt);
//      fJetCorrector->setJetA(fTR->JArea[i]);
//      fJetCorrector->setRho(fTR->Rho);
//      double correction = fJetCorrector->getCorrection();
//      jpt*=correction;
//      jenergy*=correction;
//      jpx*=correction;
//      jpy*=correction;
//      jpz*=correction;
//      nEvent.CorrectionRatio[nEvent.pfJetNum]=correction/jesC;
//      jesC=correction;

      TLorentzVector aJet(jpx,jpy,jpz,jenergy);
      
      // lepton-jet cleaning
      if ( fFullCleaning_ ) { 
        // Remove jet close to any lepton
        bool isClean(true);
        for ( size_t ilep = 0; ilep<sortedGoodLeptons.size(); ++ilep )
          if ( aJet.DeltaR(sortedGoodLeptons[ilep].p)<DRmax) isClean=false;
        if ( !isClean ) continue;
        counters[PJ].fill("... pass full lepton cleaning");
      } else {
        // Remove jet close to leptons from Z candidate
        if(aJet.DeltaR(sortedGoodLeptons[PosLepton1].p)<DRmax) continue;
        counters[PJ].fill("... pass lepton 1 veto");
        if(aJet.DeltaR(sortedGoodLeptons[PosLepton2].p)<DRmax) continue;
        counters[PJ].fill("... pass lepton 2 veto");
      }
      
      //Get Uncertainty
      fJECUnc->setJetEta(jeta);
      fJECUnc->setJetPt(jpt); // here you must use the CORRECTED jet pt
      
      float unc = fJECUnc->getUncertainty(true); 
      //removing jet from residual met
      ResidualMet=ResidualMet+aJet;
      Type1ResidualMet=Type1ResidualMet+aJet;
      //adding jet to varied MET
      float factor=(jesC+unc)/jesC;
      VariedMet=VariedMet-factor*aJet;
      Type1VariedMet=Type1VariedMet-factor*aJet;

      // Keep jets over min. pt threshold
      if ( !(jpt>20) ) continue;
      counters[PJ].fill("... pt>20.");
      
      
      if(fabs(jeta)>3.0 && fabs(jeta)<5.0) {
	if(jpt>30) nEvent.pfJetGoodNum30Fwd++;
	if(jpt>40) nEvent.pfJetGoodNum40Fwd++;
	if(jpt>50) nEvent.pfJetGoodNum50Fwd++;
	if(jpt>60) nEvent.pfJetGoodNum60Fwd++;
      }
      
      // Keep central jets
      if ( !(fabs(jeta)<3.0 ) ) continue;
      counters[PJ].fill("... |eta|<3.0");
      
      // Flag good jets failing ID
      if (!isJetID) { 
        nEvent.badJet = 1;
      } else {
        counters[PJ].fill("... pass Jet ID");
      }
      
      lepton tmpLepton;
      tmpLepton.p = aJet;
      tmpLepton.charge = 0;
      tmpLepton.index = i;
      tmpLepton.type = -1;
        
      if(!nEvent.badJet) counters[PJ].fill("... pass Jet ID");

      nEvent.pfJetPt[nEvent.pfJetNum]    = jpt;
      nEvent.pfJetEta[nEvent.pfJetNum]   = jeta;
      nEvent.pfJetPhi[nEvent.pfJetNum]   = jphi;
      nEvent.pfJetScale[nEvent.pfJetNum] = jesC;
      nEvent.pfJetScaleUnc[nEvent.pfJetNum] = unc;
      nEvent.pfJetID[nEvent.pfJetNum]    = isJetID;
      nEvent.pfJetDphiMet[nEvent.pfJetNum] = aJet.DeltaPhi(pfMETvector);
      nEvent.pfJetDphiZ[nEvent.pfJetNum] = aJet.DeltaPhi(zVector);
      nEvent.pfJetNum = nEvent.pfJetNum +1;
      nEvent.pfHT    += jpt;
      nEvent.pfGoodHT += jpt;
      sumOfPFJets += aJet;
      pfGoodJets.push_back(tmpLepton);

      if ( jpt>40 ) {
        counters[PJ].fill("... pass tight jet selection");
        nEvent.pfTightHT += jpt;
        nEvent.pfJetGoodPt[nEvent.pfJetGoodNum40]  = jpt;
        nEvent.pfJetGoodEta[nEvent.pfJetGoodNum40] = jeta;
        nEvent.pfJetGoodPhi[nEvent.pfJetGoodNum40] = jphi;
        nEvent.pfJetGoodE[nEvent.pfJetGoodNum40] = jenergy;
        nEvent.pfJetGoodID[nEvent.pfJetGoodNum40]  = isJetID;
        nEvent.bTagProbCSVBP[nEvent.pfJetGoodNum40] = fTR->JnewPFCombinedSecondaryVertexBPFJetTags[i];
        nEvent.bTagProbCSVMVA[nEvent.pfJetGoodNum40] = fTR->JnewPFCombinedSecondaryVertexMVABPFJetTags[i];
        nEvent.pfJetGoodTracks[nEvent.pfJetGoodNum40] = ntracksch;
        nEvent.pfJetGoodTracksN[nEvent.pfJetGoodNum40] = ntracksne;
        nEvent.minc+=jpt;

        if(isJetID>0) {
	        nEvent.pfJetGoodNumID++;
	        if(nEvent.bTagProbCSVBP[nEvent.pfJetGoodNum40] > 0.679) nEvent.pfJetGoodNumIDBtag++;
	}
	if(nEvent.bTagProbCSVBP[nEvent.pfJetGoodNum40] > 0.679) { 
	   nEvent.pfJetGoodPtBtag[nEvent.pfJetGoodNumBtag40]  = jpt;
           nEvent.pfJetGoodEtaBtag[nEvent.pfJetGoodNumBtag40] = jeta;
           nEvent.pfJetGoodPhiBtag[nEvent.pfJetGoodNumBtag40] = jphi;
           nEvent.pfJetGoodEBtag[nEvent.pfJetGoodNumBtag40] = jenergy;
           nEvent.pfJetGoodIDBtag[nEvent.pfJetGoodNumBtag40]  = isJetID;
	   nEvent.pfBJetDphiZ[nEvent.pfJetGoodNumBtag40]  = aJet.DeltaPhi(zVector);
           nEvent.pfJetGoodTracksBtag[nEvent.pfJetGoodNumBtag40] = ntracksch;
           nEvent.pfJetGoodTracksNBtag[nEvent.pfJetGoodNumBtag40] = ntracksne;
	   nEvent.pfJetGoodNumBtag40++;
        }
        nEvent.pfJetGoodNum40++;
      }
      if ( jpt*(jesC+unc)/jesC>30 )  nEvent.pfJetGoodNump1sigma++;
      if ( jpt*(jesC-unc)/jesC>30 )  nEvent.pfJetGoodNumn1sigma++;
      if ( jpt*(jesC+unc)/jesC>40 )  nEvent.pfJetGoodNum40p1sigma++;
      if ( jpt*(jesC-unc)/jesC>40 )  nEvent.pfJetGoodNum40n1sigma++;
      if ( jpt*(jesC+unc)/jesC>50 )  nEvent.pfJetGoodNum50p1sigma++;
      if ( jpt*(jesC-unc)/jesC>50 )  nEvent.pfJetGoodNum50n1sigma++;

      
      if ( jpt>30. && fTR->JnewPFCombinedSecondaryVertexBPFJetTags[i] > 0.679 && isJetID && abs(jeta)<2.4)  nEvent.pfJetGoodNumBtag30++;
      if ( jpt>30. && fTR->JnewPFCombinedSecondaryVertexBPFJetTags[i] > 0.898 && isJetID && abs(jeta)<2.4)  nEvent.pfJetGoodNumBtag30_Tight++;
      if ( jpt>40. && fTR->JnewPFCombinedSecondaryVertexBPFJetTags[i] > 0.898 && isJetID && abs(jeta)<2.4)  nEvent.pfJetGoodNumBtag40_Tight++;
      if ( jpt>30. && fTR->JnewPFCombinedSecondaryVertexBPFJetTags[i] > 0.244 && isJetID && abs(jeta)<2.4)  nEvent.pfJetGoodNumBtag30_Loose++;
      if ( jpt>40. && fTR->JnewPFCombinedSecondaryVertexBPFJetTags[i] > 0.244 && isJetID && abs(jeta)<2.4)  nEvent.pfJetGoodNumBtag40_Loose++;
      
      if ( jpt>30. )  nEvent.pfJetGoodNum30++;
      if ( jpt>50. )  nEvent.pfJetGoodNum50++;
      if ( jpt>60. )  nEvent.pfJetGoodNum60++;
    }
    
    
    
    if(nEvent.ZbCHS3010_pfJetGoodNum>0) {
      nEvent.ZbCHS3010_alpha=nEvent.ZbCHS3010_pfJetGoodPt[1]/nEvent.pt;
      nEvent.ZbCHS3010_alphaL5=nEvent.ZbCHS3010_L5corr*nEvent.ZbCHS3010_pfJetGoodPt[1]/nEvent.pt;
      nEvent.ZbCHS3010_alphaUp=nEvent.ZbCHS3010_alphaUp/nEvent.pt;
      nEvent.ZbCHS3010_alphaDown=nEvent.ZbCHS3010_alphaDown/nEvent.pt;
    }
    if(nEvent.ZbCHS1010_pfJetGoodNum>0) {
      nEvent.ZbCHS1010_alpha=nEvent.ZbCHS1010_pfJetGoodPt[1]/nEvent.pt;
      nEvent.ZbCHS1010_alphaUp=nEvent.ZbCHS1010_alphaUp/nEvent.pt;
      nEvent.ZbCHS1010_alphaDown=nEvent.ZbCHS1010_alphaDown/nEvent.pt;
    }
    nEvent.mpf=1+(Cleantype1METvector.Vect()*zVector.Vect())/(zVector.Px()*zVector.Px()+zVector.Py()*zVector.Py());
    nEvent.fake_mpf=1+((Cleantype1METvector.Vect()+( 1.0 / 1.05  -  1 )*zVector.Vect())*zVector.Vect())/(1.05*zVector.Px()*zVector.Px()+zVector.Py()*zVector.Py());//we simulate a missing correction of 10%
    
    
    TLorentzVector leadingJet(0,0,0,0);
    leadingJet.SetPtEtaPhiE(nEvent.pfJetGoodPt[0], nEvent.pfJetGoodEta[0], nEvent.pfJetGoodPhi[0], nEvent.pfJetGoodE[0]);
    for(int jcounter = 0; jcounter < nEvent.pfJetGoodNum40; ++jcounter) {
      TLorentzVector j1(0,0,0,0);
      j1.SetPtEtaPhiE(nEvent.pfJetGoodPt[jcounter], nEvent.pfJetGoodEta[jcounter], nEvent.pfJetGoodPhi[jcounter], nEvent.pfJetGoodE[jcounter]);
      nEvent.pfJetGoodMl[jcounter] = (sortedGoodLeptons[PosLepton1].p+sortedGoodLeptons[PosLepton2].p+j1).M();
      nEvent.pfJetGoodPtl[jcounter] = (sortedGoodLeptons[PosLepton1].p+sortedGoodLeptons[PosLepton2].p+j1).Pt();
      nEvent.pfJetGoodM[jcounter] = (j1+leadingJet).M();
    
    }

    TLorentzVector leadingBJet(0,0,0,0);
    leadingBJet.SetPtEtaPhiE(nEvent.pfJetGoodPtBtag[0], nEvent.pfJetGoodEtaBtag[0], nEvent.pfJetGoodPhiBtag[0], nEvent.pfJetGoodEBtag[0]);
    for(int jcounter = 0; jcounter < nEvent.pfJetGoodNumBtag40; ++jcounter) {
      TLorentzVector j1(0,0,0,0);
      j1.SetPtEtaPhiE(nEvent.pfJetGoodPtBtag[jcounter], nEvent.pfJetGoodEtaBtag[jcounter], nEvent.pfJetGoodPhiBtag[jcounter], nEvent.pfJetGoodEBtag[jcounter]);
      nEvent.pfJetGoodMlBtag[jcounter] = (sortedGoodLeptons[PosLepton1].p+sortedGoodLeptons[PosLepton2].p+j1).M();
      nEvent.pfJetGoodPtlBtag[jcounter] = (sortedGoodLeptons[PosLepton1].p+sortedGoodLeptons[PosLepton2].p+j1).Pt();
      nEvent.pfJetGoodMBtag[jcounter] = (j1+leadingBJet).M();

    }

    
  for(int jl=0;jl<sortedGoodLeptons.size();jl++) {
    float factor=1.0;
    if(sortedGoodLeptons[jl].type==0) {
      //electrons
      if(abs(sortedGoodLeptons[jl].p.Eta())<1.479) factor=1.006; //barrel
      else factor=1.015; // endcap
    } else {
      //muons
      factor=1.002;
    }
    
    //removing lepton from residual met
    ResidualMet=ResidualMet+sortedGoodLeptons[jl].p;
    Type1ResidualMet=Type1ResidualMet+sortedGoodLeptons[jl].p;
    //adding jet to varied MET
    VariedMet=VariedMet-factor*sortedGoodLeptons[jl].p;
    Type1VariedMet=Type1VariedMet-factor*sortedGoodLeptons[jl].p;
  }
  
  ResidualMet=1.1*ResidualMet;//all the remaining stuff is assumed to be unclustered energy
  Type1ResidualMet=1.1*Type1ResidualMet;
  
  //adding rest back in
  VariedMet+=ResidualMet;
  Type1VariedMet+=Type1ResidualMet;

  TLorentzVector s1 = sortedGoodLeptons[PosLepton1].p;
  TLorentzVector s2 = sortedGoodLeptons[PosLepton2].p;

  nEvent.met[RAW]=fTR->RawMET;
  nEvent.met[T1PFMET]=fTR->PFType1MET;
  nEvent.met[TCMET]=fTR->TCMET;
  nEvent.met[MUJESCORRMET]=fTR->MuJESCorrMET;
  nEvent.met[PFMET]=fTR->PFMET;
  nEvent.met[SUMET]=fTR->SumEt;
  
  nEvent.met[4] = nEvent.fact * nEvent.met[4];
  nEvent.pfMET  = nEvent.met[4];
  nEvent.t1pfMET= nEvent.fact * fTR->PFType1MET;

  nEvent.metUncertainty=abs((VariedMet.Pt()-nEvent.met[4])/nEvent.met[4]);
  nEvent.type1metUncertainty=abs((Type1VariedMet.Pt()-nEvent.met[1])/nEvent.met[1]);
  
  TLorentzVector pfJetVector(0,0,0,0); // for constructing SumJPt from pf jets, as Pablo
  TLorentzVector pfNoCutsJetVector(0,0,0,0); // for constructing SumJPt from pfmet (unclustered)
  TLorentzVector type1NoCutsJetVector(0,0,0,0); // same as pf, but type1 corrected
  TLorentzVector tcNoCutsJetVector(0,0,0,0); // for constructing SumJPt from tcmet (unclustered), new
  nEvent.metPhi[RAW]=0.;//kicked! caloMETvector.Phi();
  nEvent.metPhi[T1PFMET]=type1METvector.Phi();
  nEvent.metPhi[TCMET]=tcMETvector.Phi();
  nEvent.metPhi[MUJESCORRMET]=0.;
  nEvent.metPhi[PFMET]=pfMETvector.Phi();
  nEvent.metPhi[SUMET]=0.;

  // remove the leptons from PFMET and tcMET
  pfNoCutsJetVector = -pfMETvector - s1 - s2;
  type1NoCutsJetVector = -type1METvector - s1 - s2;
  tcNoCutsJetVector = -tcMETvector - s1 - s2;

  // #--- different versions of JZB
  nEvent.dphi_sumJetVSZ[TYPEONECORRPFMETJZB] = type1NoCutsJetVector.DeltaPhi(s1+s2);
  nEvent.sumJetPt[TYPEONECORRPFMETJZB] = type1NoCutsJetVector.Pt();
  nEvent.jzb[TYPEONECORRPFMETJZB] = type1NoCutsJetVector.Pt() - (s1+s2).Pt();
  nEvent.sjzb[TYPEONECORRPFMETJZB] = GausRandom(nEvent.jzb[TYPEONECORRPFMETJZB]+3.9,7); // to be used with pfMET
    
  nEvent.dphi_sumJetVSZ[PFJZB] = pfNoCutsJetVector.DeltaPhi(s1+s2); 
  nEvent.sumJetPt[PFJZB] = pfNoCutsJetVector.Pt(); 
  nEvent.jzb[PFJZB] = pfNoCutsJetVector.Pt() - (s1+s2).Pt(); // to be used with pfMET
  nEvent.sjzb[PFJZB] = GausRandom(nEvent.jzb[PFJZB]+3.9,7); // to be used with pfMET

  nEvent.dphi_sumJetVSZ[RECOILJZB] = 0.; // kicked recoil.DeltaPhi(s1+s2);
  nEvent.sumJetPt[RECOILJZB] = 0.;//kicked recoil.Pt(); 
  nEvent.jzb[RECOILJZB] = 0.;//kicked recoil.Pt() - (s1+s2).Pt(); // to be used recoil met (recoilpt[0])    
  nEvent.jzb[PFRECOILJZB] = sumOfPFJets.Pt() - (s1+s2).Pt(); // to be used recoil met (recoilpt[0])
  nEvent.sumJetPt[PFRECOILJZB] = sumOfPFJets.Pt();

  nEvent.dphi_sumJetVSZ[TCJZB] = tcNoCutsJetVector.DeltaPhi(s1+s2); // tcJZB
  nEvent.sumJetPt[TCJZB] = tcNoCutsJetVector.Pt(); 
  nEvent.jzb[TCJZB] = tcNoCutsJetVector.Pt() - (s1+s2).Pt(); // to be used with tcMET

  // --- recoil met and pf recoil met
  nEvent.met[PFRECOILMET] = (sumOfPFJets + s1 + s2).Pt(); 
  nEvent.met[RECOILMET] = 0.;//kicked (recoil + s1 + s2).Pt();
    
  // Statistics ///////////////////////////////////////
  string type("");
  switch ( (nEvent.id1+1)*(nEvent.id2+1) ) {
  case 1: type = "ee"; break;
  case 2: type = "em"; break;
  case 4: type = "mm"; break;
  default: type = "unknown";
  }
  counters[EV].fill("... "+type+" pairs");     
  if ( nEvent.pfJetGoodNum40>= 2 ) {
    counters[EV].fill("... "+type+" + 2 jets");
    if ( fabs(nEvent.mll-91)<20 ) {
      counters[EV].fill("... "+type+" + 2 jets + require Z");
      if ( nEvent.jzb[1]>50 ) {
        counters[EV].fill("... "+type+" + 2 jets + require Z + JZB>50");
      }
    }
  }
  // Trigger information
  map<string,int>::iterator itend = fHLTLabelMap.end();
  char buf[256];
  counters[TR].fill("All selected events");
  for ( map<string,int>::iterator it = fHLTLabelMap.begin(); it != itend; ++it ) {
    int bit = it->second;
    bool passed = fTR->HLTResults[bit];
    if ( passed ) {
      //sprintf(buf,"... %s (%02d)",(it->first).c_str(),fTR->HLTPrescale[bit]);
      counters[TR].fill( (it->first), fTR->HLTResults[bit] );
    }
  }
  ////////////////////////////////////////////////////



  // --- store number of good leptons in the event 
  nEvent.leptonNum = int(sortedGoodLeptons.size());
  for ( size_t i=0; i<sortedGoodLeptons.size(); ++i ) {
    TLorentzVector lp(sortedGoodLeptons[i].p);
    nEvent.leptonPt[i] = lp.Pt();
    nEvent.leptonEta[i] = lp.Eta();
    nEvent.leptonPhi[i] = lp.Phi();
    nEvent.leptonCharge[i] = sortedGoodLeptons[i].charge;
    nEvent.leptonId[i] = sortedGoodLeptons[i].type ;
      

    for(size_t j=i+1; j<sortedGoodLeptons.size();j++) // store lepton pair masses
      {
        TLorentzVector lp1(sortedGoodLeptons[i].p);
        TLorentzVector lp2(sortedGoodLeptons[j].p);
        int old_id1 = (sortedGoodLeptons[i].type+1)*sortedGoodLeptons[i].charge;
        int old_id2 = (sortedGoodLeptons[j].type+1)*sortedGoodLeptons[j].charge;
        if(nEvent.leptonPairNum<jMax)
          {
            nEvent.leptonPairMass[nEvent.leptonPairNum] = (lp1+lp2).M();
            nEvent.leptonPairDphi[nEvent.leptonPairNum] = lp1.DeltaPhi(lp2);
            nEvent.leptonPairId[nEvent.leptonPairNum] = old_id1*old_id2;
            nEvent.leptonPairNum=nEvent.leptonPairNum+1;
          }
      }
  }

  if(nEvent.leptonNum>=3 && TriLepton1<sortedGoodLeptons.size() && TriLepton2<sortedGoodLeptons.size() && TriLepton3<sortedGoodLeptons.size() ) {
    TLorentzVector lp1(sortedGoodLeptons[TriLepton1].p);
    TLorentzVector lp2(sortedGoodLeptons[TriLepton2].p);
    TLorentzVector lp3(sortedGoodLeptons[TriLepton3].p);
    
    nEvent.tri_pt1=sortedGoodLeptons[TriLepton1].p.Pt();
    nEvent.tri_pt2=sortedGoodLeptons[TriLepton2].p.Pt();
    nEvent.tri_pt3=sortedGoodLeptons[TriLepton3].p.Pt();
    
    nEvent.tri_id1=sortedGoodLeptons[TriLepton1].type;
    nEvent.tri_id2=sortedGoodLeptons[TriLepton2].type;
    nEvent.tri_id3=sortedGoodLeptons[TriLepton3].type;
    
    nEvent.tri_ch1=sortedGoodLeptons[TriLepton1].charge;
    nEvent.tri_ch2=sortedGoodLeptons[TriLepton2].charge;
    nEvent.tri_ch3=sortedGoodLeptons[TriLepton3].charge;
    
    nEvent.tri_mlll=(lp1+lp2+lp3).M();
    nEvent.tri_mll=(lp1+lp2).M();
    nEvent.tri_submll=(lp2+lp3).M();
    
    nEvent.tri_badmll=(sortedGoodLeptons[BadTriLepton1].p+sortedGoodLeptons[BadTriLepton2].p).M();
    nEvent.tri_badsubmll=(sortedGoodLeptons[BadTriLepton2].p+sortedGoodLeptons[BadTriLepton3].p).M();
    
    float AngleBetweenMETandThirdLepton=sortedGoodLeptons[TriLepton3].p.DeltaPhi(pfMETvector);
    nEvent.tri_mT = sqrt(2 * nEvent.met[4] * sortedGoodLeptons[TriLepton3].p.Pt() * ( 1 - cos(AngleBetweenMETandThirdLepton)));
    nEvent.tri_index1=TriLepton1;
    nEvent.tri_index2=TriLepton2;
    nEvent.tri_index3=TriLepton3;
    
    AngleBetweenMETandThirdLepton=sortedGoodLeptons[BadTriLepton3].p.DeltaPhi(pfMETvector);
    nEvent.tri_badmT = sqrt(2 * nEvent.met[4] * sortedGoodLeptons[BadTriLepton3].p.Pt() * ( 1 - cos(AngleBetweenMETandThirdLepton)));
    nEvent.tri_badindex1=BadTriLepton1;
    nEvent.tri_badindex2=BadTriLepton2;
    nEvent.tri_badindex3=BadTriLepton3;

    nEvent.tri_badid1=sortedGoodLeptons[BadTriLepton1].type;
    nEvent.tri_badid2=sortedGoodLeptons[BadTriLepton2].type;
    nEvent.tri_badid3=sortedGoodLeptons[BadTriLepton3].type;
    
    int i1 = sortedGoodLeptons[TriLepton1].index;
    int i2 = sortedGoodLeptons[TriLepton2].index;
    int i3 = sortedGoodLeptons[TriLepton3].index;

    int bi1 = sortedGoodLeptons[BadTriLepton1].index;
    int bi2 = sortedGoodLeptons[BadTriLepton2].index;
    int bi3 = sortedGoodLeptons[BadTriLepton3].index;

    
    if(fDataType_ == "mc" ) {
      nEvent.tri_genMID1 = (sortedGoodLeptons[TriLepton1].type?fTR->MuGenMID[i1]:fTR->ElGenMID[i1]);
      nEvent.tri_genMID2 = (sortedGoodLeptons[TriLepton2].type?fTR->MuGenMID[i2]:fTR->ElGenMID[i2]);
      nEvent.tri_genMID3 = (sortedGoodLeptons[TriLepton3].type?fTR->MuGenMID[i3]:fTR->ElGenMID[i3]);
      
      if(abs(nEvent.tri_genMID1)==15) nEvent.tri_genMID1=(sortedGoodLeptons[TriLepton1].type?fTR->MuGenGMID[i1]:fTR->ElGenGMID[i1]);
      if(abs(nEvent.tri_genMID2)==15) nEvent.tri_genMID2=(sortedGoodLeptons[TriLepton2].type?fTR->MuGenGMID[i2]:fTR->ElGenGMID[i2]);
      if(abs(nEvent.tri_genMID3)==15) nEvent.tri_genMID3=(sortedGoodLeptons[TriLepton3].type?fTR->MuGenGMID[i3]:fTR->ElGenGMID[i3]);
      
      nEvent.tri_badgenMID1 = (sortedGoodLeptons[BadTriLepton1].type?fTR->MuGenMID[bi1]:fTR->ElGenMID[bi1]);
      nEvent.tri_badgenMID2 = (sortedGoodLeptons[BadTriLepton2].type?fTR->MuGenMID[bi2]:fTR->ElGenMID[bi2]);
      nEvent.tri_badgenMID3 = (sortedGoodLeptons[BadTriLepton3].type?fTR->MuGenMID[bi3]:fTR->ElGenMID[bi3]);
      
      if(abs(nEvent.tri_badgenMID1)==15) nEvent.tri_badgenMID1=(sortedGoodLeptons[BadTriLepton1].type?fTR->MuGenGMID[bi1]:fTR->ElGenGMID[bi1]);
      if(abs(nEvent.tri_badgenMID2)==15) nEvent.tri_badgenMID2=(sortedGoodLeptons[BadTriLepton2].type?fTR->MuGenGMID[bi2]:fTR->ElGenGMID[bi2]);
      if(abs(nEvent.tri_badgenMID3)==15) nEvent.tri_badgenMID3=(sortedGoodLeptons[BadTriLepton3].type?fTR->MuGenGMID[bi3]:fTR->ElGenGMID[bi3]);
      
      nEvent.tri_GoodZMatch=false;
      nEvent.tri_GoodWMatch=false;
      if(nEvent.tri_genMID1==23 && nEvent.tri_genMID1 ==23) nEvent.tri_GoodZMatch=true;
      if(abs(nEvent.tri_genMID3)==24 ) nEvent.tri_GoodWMatch=true;
    }

    nEvent.tri_dR12 = lp1.DeltaR(lp2);
    nEvent.tri_dR13 = lp1.DeltaR(lp3); 
    nEvent.tri_dR23 = lp2.DeltaR(lp3);

    nEvent.tri_eta1 = lp1.Eta();
    nEvent.tri_eta2 = lp2.Eta();
    nEvent.tri_eta3 = lp3.Eta();
    
  }
  nEvent.dphiZpfMet = (s1+s2).DeltaPhi(pfMETvector);
  nEvent.dphiZs1 = (s1+s2).DeltaPhi(s1);
  nEvent.dphiZs2 = (s1+s2).DeltaPhi(s2);
  nEvent.dphiMet1 = sortedGoodLeptons[PosLepton1].p.DeltaPhi(pfMETvector);
  nEvent.dphiMet2 = sortedGoodLeptons[PosLepton2].p.DeltaPhi(pfMETvector);
  nEvent.dphitcMet1 = sortedGoodLeptons[PosLepton1].p.DeltaPhi(tcNoCutsJetVector);
  nEvent.dphitcMet2 = sortedGoodLeptons[PosLepton2].p.DeltaPhi(tcNoCutsJetVector);
  nEvent.dphipft1Met1 = sortedGoodLeptons[PosLepton1].p.DeltaPhi(type1NoCutsJetVector);
  nEvent.dphipft1Met2 = sortedGoodLeptons[PosLepton2].p.DeltaPhi(type1NoCutsJetVector);
  nEvent.dphipfRecoilMet1 = sortedGoodLeptons[PosLepton1].p.DeltaPhi(-sumOfPFJets - s1 - s2); // pf recoil met
  nEvent.dphipfRecoilMet2 = sortedGoodLeptons[PosLepton2].p.DeltaPhi(-sumOfPFJets - s1 - s2); // pf recoil met
    
  // Store minimum dphi between some mets and any kind of lepton
  for ( size_t i=0; i<sortedGoodLeptons.size(); ++i ) {
    TLorentzVector lp(sortedGoodLeptons[i].p);
    if ( fabs(pfMETvector.DeltaPhi(lp))<fabs(nEvent.dphiMetLep[PFMET]) ) nEvent.dphiMetLep[PFMET] = pfMETvector.DeltaPhi(lp);
    if ( fabs((sumOfPFJets + s1 + s2).DeltaPhi(lp))< fabs(nEvent.dphiMetLep[PFRECOILMET]) ) nEvent.dphiMetLep[PFRECOILMET] = (sumOfPFJets + s1 + s2).DeltaPhi(lp);
  }

  // Store minimum dphi between some mets and any good jet
  for ( size_t i=0; i<pfGoodJets.size(); ++i ) {
    TLorentzVector jp(pfGoodJets[i].p);
    if ( fabs(pfMETvector.DeltaPhi(jp))<fabs(nEvent.dphiMetJet[PFMET]) )
      nEvent.dphiMetJet[PFMET] = pfMETvector.DeltaPhi(jp);
  }
  nEvent.dphiMetSumJetPt[PFMET] = pfNoCutsJetVector.DeltaPhi(pfMETvector);

  // Store some additional MET information
  nEvent.metPar[PFMET]  = pfMETvector.Dot(s1+s2);
  nEvent.metPerp[PFMET] = pfMETvector.Perp((s1+s2).Vect());
    
  // Store some generator information on selected leptons
  if ( isMC ) {
    TLorentzVector GenMETvector(fTR->GenMETpx,fTR->GenMETpy,0,0);
    int i1 = sortedGoodLeptons[PosLepton1].index;
    int i2 = sortedGoodLeptons[PosLepton2].index;

    TLorentzVector genLep1; 
    if ( sortedGoodLeptons[PosLepton1].type )
      genLep1.SetPtEtaPhiE(fTR->MuGenPt[i1],fTR->MuGenEta[i1],fTR->MuGenPhi[i1],fTR->MuGenE[i1]);
    else
      genLep1.SetPtEtaPhiE(fTR->ElGenPt[i1],fTR->ElGenEta[i1],fTR->ElGenPhi[i1],fTR->ElGenE[i1]);
    TLorentzVector genLep2;
    if ( sortedGoodLeptons[PosLepton2].type )
      genLep2.SetPtEtaPhiE(fTR->MuGenPt[i2],fTR->MuGenEta[i2],fTR->MuGenPhi[i2],fTR->MuGenE[i2]);
    else
      genLep2.SetPtEtaPhiE(fTR->ElGenPt[i2],fTR->ElGenEta[i2],fTR->ElGenPhi[i2],fTR->ElGenE[i2]);
      
    nEvent.genRecoilSel = (-GenMETvector - genLep1 - genLep2).Pt();
    nEvent.genZPtSel    = (genLep1 + genLep2).Pt();
    nEvent.genMllSel    = (genLep1 + genLep2).M();

    nEvent.genMID1     = (sortedGoodLeptons[PosLepton1].type?fTR->MuGenMID[i1]:fTR->ElGenMID[i1]); // WW study
    nEvent.genMID2     = (sortedGoodLeptons[PosLepton2].type?fTR->MuGenMID[i2]:fTR->ElGenMID[i2]); // WW study

    nEvent.genGMID1    = (sortedGoodLeptons[PosLepton1].type?fTR->MuGenGMID[i1]:fTR->ElGenGMID[i1]); // WW study
    nEvent.genGMID2    = (sortedGoodLeptons[PosLepton2].type?fTR->MuGenGMID[i2]:fTR->ElGenGMID[i2]); // WW study

    nEvent.genPt1Sel    = genLep1.Pt();
    nEvent.genPt2Sel    = genLep2.Pt();
    nEvent.genEta1Sel   = genLep1.Eta();
    nEvent.genEta2Sel   = genLep2.Eta();
    nEvent.genId1Sel    = (sortedGoodLeptons[PosLepton1].type?fTR->MuGenID[i1]:fTR->ElGenID[i1]);
    nEvent.genId2Sel    = (sortedGoodLeptons[PosLepton2].type?fTR->MuGenID[i2]:fTR->ElGenID[i2]);
    nEvent.genJZBSel    = nEvent.genRecoilSel - (genLep1 + genLep2).Pt();
  }
  myTree->Fill();
}

void JZBAnalysis::End(TFile *f){
  f->cd();	

  myTree->Write();
  weight_histo->Write();

  // Dump statistics
  if (1) { // Put that to 0 if you are annoyed
    std::cout << setfill('=') << std::setw(70) << "" << std::endl;
    std::cout << "Statistics" << std::endl;
    std::cout << setfill('-') << std::setw(70) << "" << setfill(' ') << std::endl;
    for ( counters_t iCount=count_begin; iCount<count_end; 
          iCount = counters_t(iCount+1) ) {
      counters[iCount].print();
    }
  }
  
}

template<class T>
std::string JZBAnalysis::any2string(T i)
{
  std::ostringstream buffer;
  buffer << i;
  return buffer.str();
}

const bool JZBAnalysis::IsSoftMuon(const int index) {

  if ( !fTR->MuIsTrackerMuon[index] )       return false;
  if ( !fTR->MuIsTMLSTight[index])          return false;
  if ( !(fabs(fTR->MuEta[index])<2.4) )     return false;
  if ( !(fTR->MuNSiLayers[index] > 5) )     return false;
  if ( !(fabs(fTR->MuD0PV[index]) < 0.2) )  return false;
  if ( !(fabs(fTR->MuDzPV[index]) < 0.1 ) ) return false;

  double Iso = MuPFIso(index);
  if (!(fTR->MuPt[index] <= 20 || (fTR->MuPt[index] > 20 && Iso > 0.15)) ) return false;

  return true;

}

const bool JZBAnalysis::IsCustomMu2012(const int index){

  // Basic muon cleaning and ID

  // Acceptance cuts
  if (!(fTR->MuPt[index] > 10.0) )       return false;
  counters[MU].fill(" ... pt > 10");
  if (!(fabs(fTR->MuEta[index])<2.4) ) return false;
  counters[MU].fill(" ... |eta| < 2.4");


  // Quality cuts
  if ( !fTR->MuIsGlobalMuon[index] )  return false;
  counters[MU].fill(" ... is global muon");
  if ( !fTR->MuIsTrackerMuon[index] ) return false;
  counters[MU].fill(" ... is tracker muon");
  if ( !fTR->MuIsPFMuon[index] )        return false;
  counters[MU].fill(" ... is pf muon");

  // Hits
  if ( !(fTR->MuNChi2[index] < 10) )     return false;
  counters[MU].fill(" ... nChi2 < 10");
  if ( !(fTR->MuNMuHits[index] > 0) )     return false;
  counters[MU].fill(" ... nValidHits > 0");
  if ( !(fTR->MuNPxHits[index] > 0) )       return false;
  counters[MU].fill(" ... nPxHits > 0");
  //FR if ( !(fTR->MuNMatchedStations[index] > 1) )      return false;
  if ( !(fTR->MuNMatches[index] > 1) )      return false;
  counters[MU].fill(" ... nMatches > 1");
  if ( !(fTR->MuNSiLayers[index] > 5) )      return false;
  counters[MU].fill(" ... nLayers > 5");


  // Vertex compatibility
  if ( !(fabs(fTR->MuD0PV[index]) < 0.02) ) return false; //still open
  counters[MU].fill(" ... D0(pv) < 0.02");
  //HPA recommendation not POG
  //FR if ( !(fabs(fTR->MuDzPV[index]) < 0.1 ) ) return false; //still open
  //FR counters[MU].fill(" ... DZ(pv) < 0.1");
  if ( !(fabs(fTR->MuDzPV[index]) < 0.2 ) ) return false; //still open
  counters[MU].fill(" ... DZ(pv) < 0.2");


  // Flat isolation below 20 GeV (only for synch.: we cut at 20...)
  double Iso = MuPFIso(index);
  if ( !(Iso < 0.15) ) return false;
  counters[MU].fill(" ... Iso < 0.15");


  return true;
}



const float JZBAnalysis::GetLeptonWeight(int id1, float pt1, float eta1, int id2, float pt2, float eta2, float &EffErr) {

  // FIXME: Need to update!
  //EffErr = 0.0;
  //return 1.0;

  float weight1;
  float error1;
  float weight2;
  float error2; 

  if(id1==0) weight1 = GetElectronWeight(eta1, pt1, error1);
  if(id2==0) weight2 = GetElectronWeight(eta2, pt2, error2);
  if(id1==1) weight1 = GetMuonWeight(eta1, pt1, error1);
  if(id2==1) weight2 = GetMuonWeight(eta2, pt2, error2);

  EffErr = TMath::Sqrt(error1*error1*weight1*weight1+error2*error2*weight2*weight2);
  return weight1*weight2;

  
  /*if(id1==id2&&id1==0) {
      EffErr=0.01;
      return 0.99;
  }
  if(id1==id2&&id1==1) {
     //mm
      EffErr=0.02;
      return 0.95;
    }
  if(id1!=id2) {
      //em
      EffErr=0.03;
      return 0.98;
   }*/

}



const float JZBAnalysis::GetMuonWeight(float eta1, float pt1, float &error1) {

  float weight1 = 1.0;

  if(abs(eta1)<1.2) {
    if(pt1 > 10 && pt1 < 20) {
      weight1 = 1.06;
      error1 = 0.01;
    } 
    else if(pt1 > 20 && pt1 < 30) {
      weight1 = 1.078;
      error1 = 0.003;
    }
    else if(pt1 > 30 && pt1 < 40) {
      weight1 = 1.041;
      error1 = 0.002;
    }
    else if(pt1 > 40 && pt1 < 60) {
      weight1 = 1.0184;
      error1 = 0.0008;
    }
    else if(pt1 > 60) { //or only between 60 and 100????
      weight1 = 1.007;
      error1 = 0.002;
    }
  }  
  else if(abs(eta1)>1.2) {
   if(pt1 > 10 && pt1 < 20) {
      weight1 = 0.970;
      error1 = 0.051;
    }
    else if(pt1 > 20 && pt1 < 30) {
      weight1 = 1.025;
      error1 = 0.003;
    }
    else if(pt1 > 30 && pt1 < 40) {
      weight1 = 1.017;
      error1 = 0.001;
    }
    else if(pt1 > 40 && pt1 < 60) {
      weight1 = 1.009;
      error1 = 0.0007;
    }
    else if(pt1 > 60) { //or only between 60 and 100???
      weight1 = 1.005;
      error1 = 0.002;
    }
  } 

  return weight1;

}



const float JZBAnalysis::GetElectronWeight(float eta1, float pt1, float &error1) {


   //This is for medium selection
   float weight1 = 1.0;

   if(abs(eta1)<0.8) {
     if(pt1 > 10.0 && pt1 < 15.0) {
       weight1 = 0.879; 
       error1 = 0.040;
     }
     else if(pt1 > 15.0 && pt1 < 20.0) {
       weight1 = 0.946; 
       error1 = 0.018;
     }
     else if(pt1 > 20.0 && pt1 < 30.0) {
       weight1 = 1.017; 
       error1 = 0.004;
     }
     else if(pt1 > 30.0 && pt1 < 40.0) {
       weight1 = 1.019; 
       error1 = 0.002;
     }
     else if(pt1 > 40.0 && pt1 < 50.0) {
       weight1 = 1.015; 
       error1 = 0.001;
     }
     else if(pt1 > 50.0) {
       weight1 = 1.005; 
       error1 = 0.002;
     }
   }
   else if(abs(eta1) > 0.8 && abs(eta1) < 1.442) {
     if(pt1 > 10.0 && pt1 < 15.0) {
       weight1 = 0.885;
       error1 = 0.043;
     }
     else if(pt1 > 15.0 && pt1 < 20.0) {
       weight1 = 0.932;
       error1 = 0.021;
     }
     else if(pt1 > 20.0 && pt1 < 30.0) {
       weight1 = 0.991;
       error1 = 0.013;
     }
     else if(pt1 > 30.0 && pt1 < 40.0) {
       weight1 = 1.002;
       error1 = 0.002;
     }
     else if(pt1 > 40.0 && pt1 < 50.0) {
       weight1 = 1.00;
       error1 = 0.002;
     }
     else if(pt1 > 50.0) {
       weight1 = 0.992;
       error1 = 0.003;
     }
   }
   else if(abs(eta1) > 1.442 && abs(eta1) < 1.556) {
     if(pt1 > 10.0 && pt1 < 15.0) {
       weight1 = 1.020;
       error1 = 0.170;
     }
     else if(pt1 > 15.0 && pt1 < 20.0) {
       weight1 = 1.045;
       error1 = 0.082;
     }
     else if(pt1 > 20.0 && pt1 < 30. ) {
       weight1 = 1.176;
       error1 = 0.025;
     }
     else if(pt1 > 30.0 && pt1 < 40.0) {
       weight1 = 1.038;
       error1 = 0.011;
     }
     else if(pt1 > 40.0 && pt1 < 50.0) {
       weight1 = 0.985;
       error1 = 0.009;
     }
     else if(pt1 > 50.0) {
       weight1 = 0.990;
       error1 = 0.037;
     }
   }
   else if(abs(eta1) > 1.556 && abs(eta1) < 2.0) {
     if(pt1 > 10.0 && pt1 < 15.0) {
       weight1 = 0.873;
       error1 = 0.073;
     }
     else if(pt1 > 15.0 && pt1 < 20.0) {
       weight1 = 0.958;
       error1 = 0.031;
     }
     else if(pt1 > 20.0 && pt1 < 30. ) {
       weight1 = 1.010;
       error1 = 0.009;
     }
     else if(pt1 > 30.0 && pt1 < 40.0) {
       weight1 = 1.010;
       error1 = 0.004;
     }
     else if(pt1 > 40.0 && pt1 < 50.0) {
       weight1 = 1.013;
       error1 = 0.003;
     }
     else if(pt1 > 50.0) {
       weight1 = 1.005;
       error1 = 0.005;
     }
   }
   else if(abs(eta1) > 2.0 && abs(eta1) < 2.5) {
     if(pt1 > 10.0 && pt1 < 15.0) {
       weight1 = 1.080;
       error1 = 0.075;
     }
     else if(pt1 > 15.0 && pt1 < 20.0) {
       weight1 = 1.175;
       error1 = 0.036;
     }
     else if(pt1 > 20.0 && pt1 < 30. ) {
       weight1 = 1.110;
       error1 = 0.005;
     }
     else if(pt1 > 30.0 && pt1 < 40.0) {
       weight1 = 1.074;
       error1 = 0.005;
     }
     else if(pt1 > 40.0 && pt1 < 50.0) {
       weight1 = 1.042;
       error1 = 0.004;
     }
     else if(pt1 > 50.0) {
       weight1 = 1.020;
       error1 = 0.007;
     }
   }

 return weight1;

}

const float JZBAnalysis::IndividualEffArea(float abseta, string type) {
    //from: https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedPhotonID2012
    
    if(type=="chargedhad") {
        if(abs(abseta)<1.0) return 	0.002;
        if(abs(abseta)<1.479) return 0.003;
        if(abs(abseta)<2.0) return 0.004;
        if(abs(abseta)<2.2) return 0.006;
        if(abs(abseta)<2.3) return 0.006;
        if(abs(abseta)<2.4) return 0.004;
        return 0.003;
    }
    
    if(type=="neutralhad") {
        if(abs(abseta)<1.0) return 	0.024;
        if(abs(abseta)<1.479) return 0.037;
        if(abs(abseta)<2.0) return 0.037;
        if(abs(abseta)<2.2) return 0.034;
        if(abs(abseta)<2.3) return 0.043;
        if(abs(abseta)<2.4) return 0.047;
        return 0.066;
    }
    
    if(type=="photon") {
        if(abs(abseta)<1.0) return 	0.053 ;
        if(abs(abseta)<1.479) return 0.052;
        if(abs(abseta)<2.0) return 0.037;
        if(abs(abseta)<2.2) return 0.073;
        if(abs(abseta)<2.3) return 0.107;
        if(abs(abseta)<2.4) return 0.123;
        return 0.133;
    }
    
    cout << "YOUR PROVIDED TYPE, " << type << " HAS NOT BEEN RECOGNIZED. RETURNING HUGE NEGATIVE EFFECTIVE AREA" << endl;
    return -10e8;
}
    
const float JZBAnalysis::EffArea(float abseta) {
  abseta=fabs(abseta); // making sure we're looking at |eta|
  if(abseta<1.0) return 0.10;
  if(abseta<1.479) return 0.12;
  if(abseta<2.0) return 0.085;
  if(abseta<2.2) return 0.11;
  if(abseta<2.3) return 0.12;
  if(abseta<2.4) return 0.12;
  return 0.13;
}

float JZBAnalysis::PhoPFIso(int index){
    float eta=fTR->PhoEta[index];
    float pfchargedhadiso = max(fTR->PhoNewIsoPFCharged[index]-IndividualEffArea(eta,"chargedhad")*fTR->RhoForIso,(float)0.0)/fTR->PhoPt[index];
    float pfneutralhadiso = max(fTR->PhoNewIsoPFNeutral[index]-IndividualEffArea(eta,"neutralhad")*fTR->RhoForIso,(float)0.0)/fTR->PhoPt[index];
    float pfphotoniso     = max(fTR->PhoNewIsoPFPhoton[index]-IndividualEffArea(eta,"photon")*fTR->RhoForIso,(float)0.0)/fTR->PhoPt[index];
    return pfchargedhadiso+pfneutralhadiso+pfphotoniso;
}

float JZBAnalysis::MuPFIso(int index){
   double neutral = (fTR->MuPfIsoR03NeHad[index] + fTR->MuPfIsoR03Photon[index] - 0.5*fTR->MuPfIsoR03SumPUPt[index] );
   float iso = ( fTR->MuPfIsoR03ChHad[index] + TMath::Max(0., neutral) ) / fTR->MuPt[index];
   return iso;
}

float JZBAnalysis::ElPFIso(int index){
   double neutral = fTR->ElEventelPFIsoValueNeutral03PFIdStandard[index] + fTR->ElEventelPFIsoValueGamma03PFIdStandard[index];
   double rhocorr = fTR->RhoForIso * EffArea(fTR->ElEta[index]);
   double iso = ( fTR->ElEventelPFIsoValueCharged03PFIdStandard[index] + TMath::Max(0., neutral - rhocorr) )/ fTR->ElPt[index];
   return iso;
}

const bool JZBAnalysis::IsCustomEl2012(const int index) {
  
  if(!(fabs(fTR->ElEta[index]) < 2.5) ) return false;
  counters[EL].fill(" ... |eta| < 2.5");

  if(!(fTR->ElPt[index]) > 10.0 ) return false;
  counters[EL].fill(" ... pT > 10");

  // Medium Working Point
  if ( fabs(fTR->ElEta[index]) < 1.479 ) { // Barrel
     if(!(fabs(fTR->ElDeltaEtaSuperClusterAtVtx[index])<0.007)) return false;
     if(!(fabs(fTR->ElDeltaPhiSuperClusterAtVtx[index])<0.15)) return false;
     if(!(fTR->ElSigmaIetaIeta[index]<0.01)) return false;
     if(!(fTR->ElHcalOverEcal[index]<0.12)) return false;
  } else { // Endcap
     if(!(fabs(fTR->ElDeltaEtaSuperClusterAtVtx[index])<0.009 )) return false;
     if(!(fabs(fTR->ElDeltaPhiSuperClusterAtVtx[index])<0.10 )) return false;
     if(!(fTR->ElSigmaIetaIeta[index]<0.03)) return false;
     if(!(fTR->ElHcalOverEcal[index]<0.10)) return false;
  }
  
  counters[EL].fill(" ... pass additional electron ID cuts");

  if(!(abs(fTR->ElD0PV[index])<0.02)) return false;
  counters[EL].fill(" ... D0(PV)<0.02");
  //FR if(!(abs(fTR->ElDzPV[index])<0.1)) return false;
  //FR counters[EL].fill(" ... DZ(PV)<0.1");
  if(!(abs(fTR->ElDzPV[index])<0.2)) return false;
  counters[EL].fill(" ... DZ(PV)<0.2");

//  if(!(fTR->ElPassConversionVeto[index])) return false;
  if(!(fTR->ElNumberOfMissingInnerHits[index]<=1)) return false;
  counters[EL].fill(" ... N(missing inner hits) <= 1");
  if(!fTR->ElPassConversionVeto[index]) return false;
  counters[EL].fill(" ... passed conversion rejection");
  

  float e=fTR->ElCaloEnergy[index];
  float p=fTR->ElCaloEnergy[index]/fTR->ElESuperClusterOverP[index];
  if(!(fabs(1/e-1/p)<0.05)) return false;
  counters[EL].fill(" ... |1/e-1/p|<0.05");
  
  // ECAL gap veto
  if ( fabs(fTR->ElSCEta[index]) > 1.4442 && fabs(fTR->ElSCEta[index]) < 1.566 )  return false;  
  counters[EL].fill(" ... not in ECAL gap");

//fbrem : fTElNBrems (reco::GsfElectron::fbrem())  --> no cut?
  
  float pfIso = ElPFIso(index);
  if ( fabs(fTR->ElEta[index]) < 1.479 || fTR->ElPt[index]>20.0) { // Barrel
    if ( !((pfIso  < 0.15) ) ) return false;
  } else {
    //Endcap with pt<20
    if ( !((pfIso  < 0.10) ) ) return false;
  }
  counters[EL].fill(" ... pfIso  < 0.15 (or 0.1 for endcaps with pt<20)");

  return true;
}

// Check if electron is from photon conversion
const bool JZBAnalysis::IsConvertedPhoton( const int eIndex ) {
 
  int elIDWP90 = fTR->ElIDsimpleWP90relIso[eIndex];
  if ( elIDWP90 < 4 ) return true;
  counters[EL].fill(" ... passes conversion rejection");
  return false;
 
}

const bool JZBAnalysis::IsCustomJet(const int index){
  // Basic Jet ID cuts (loose Jet ID)
  // See https://twiki.cern.ch/twiki/bin/view/CMS/JetID

//  if ( !(fTR->CAJID_n90Hits[index] > 1) ) return false;
  counters[JE].fill(" ... n90Hits > 1");
//  if ( !(fTR->CAJID_HPD[index] < 0.98)  ) return false;
  counters[JE].fill(" ... HPD < 0.98");

  if ( fabs(fTR->CAJEta[index])<3.0 ) {
    if ( !(fTR->CAJEMfrac[index] > 0.01)  ) return false;
  } else {
    if ( !(fTR->CAJEMfrac[index] > -0.9)  ) return false;
    if ( fTR->CAJPt[index] > 80 && !(fTR->CAJEMfrac[index]<1) ) return false;
  }
  counters[JE].fill(" ... pass EMfrac cut");

  return true;
}

void JZBAnalysis::GeneratorInfo(void) {
  // Try to find an Z->ll pair inside the acceptance
  double minPt = 5.;
  double mllCut = 20.;
  double maxEta = 2.4;
  double minJPt = 40;
  double maxJEta = 3.0;
  // First, look for leptons in acceptance
  vector<lepton> gLeptons;
  for ( int gIndex=0;gIndex<fTR->NGenLeptons; ++gIndex ) {
    if ( fTR->GenLeptonPt[gIndex]>minPt && 
         fabs(fTR->GenLeptonEta[gIndex])<maxEta &&
         ( abs(fTR->GenLeptonID[gIndex])==11 ||
           abs(fTR->GenLeptonID[gIndex])==13 ) 
         )
      {
        TLorentzVector tmpVector;
        tmpVector.SetPtEtaPhiM(fTR->GenLeptonPt[gIndex],fTR->GenLeptonEta[gIndex], 
                               fTR->GenLeptonPhi[gIndex],0.);
        lepton tmpLepton;
        tmpLepton.p      = tmpVector;
        tmpLepton.charge = fTR->GenLeptonID[gIndex]/abs(fTR->GenLeptonID[gIndex]);
        tmpLepton.index  = gIndex;
        tmpLepton.type   = fTR->GenLeptonID[gIndex];
        tmpLepton.genPt  = tmpVector.Pt();
        gLeptons.push_back(tmpLepton); 
        //if ( fTR->GenLeptonMID[gIndex] ==23 ) gLeptons.push_back(tmpLepton); // WW study
      }         
  }
  
  // Gen leptons are not sorted by Pt...
  vector<lepton> sortedGLeptons = sortLeptonsByPt(gLeptons);

  // Store actual number of leptons passing selection
  nEvent.genNleptons = gLeptons.size();

  // Now fill information
  TLorentzVector GenMETvector(fTR->GenMETpx,fTR->GenMETpy,0,0);
  nEvent.genMET     = fTR->GenMET;

  // Number of good jets
  nEvent.genNjets = 0;
  TLorentzVector sumOfGenJets(0,0,0,0);
  for ( int jIndex=0; jIndex<fTR->NGenJets; ++jIndex) {
    if ( fTR->GenJetPt[jIndex]<minJPt ) continue;
    if ( fabs(fTR->GenJetEta[jIndex])>maxJEta ) continue;
    ++nEvent.genNjets;
    if ( fabs(fTR->GenJetEta[jIndex])>3.0 ) continue;
    TLorentzVector gJet;
    gJet.SetPtEtaPhiE(fTR->GenJetPt[jIndex],fTR->GenJetEta[jIndex],fTR->GenJetPhi[jIndex],fTR->GenJetE[jIndex]);
    sumOfGenJets+=gJet;
  }


  //FR: Deactivate Z matching
  size_t i1 = 0, i2 = 1;

   for(; i2<sortedGLeptons.size();i2++) {
     if(sortedGLeptons[i1].charge*sortedGLeptons[i2].charge<0) break;
   }
   
//   // Select the two highest-pt leptons compatible with a Z
//   size_t i1 = 0, i2 = 0;
//   if ( sortedGLeptons.size()>1 ) {
//     for ( size_t i=0; i<sortedGLeptons.size()-1; ++i ) {
//       i1 = i;
//       TLorentzVector lp1(sortedGLeptons[i].p);
//       for ( size_t j=i+1; j<sortedGLeptons.size(); ++j ) {
//         i2 = j;
//         TLorentzVector lp2(sortedGLeptons[j].p);
//         if ( fabs( (lp1+lp2).M() - 91.2 ) < mllCut ) break;
//       }
//     }
//   }
  
  if(sortedGLeptons.size()>0)
    {
      nEvent.genPt1     = sortedGLeptons[i1].p.Pt();
      nEvent.genId1     = sortedGLeptons[i1].type;
      nEvent.genPhi1    = sortedGLeptons[i1].p.Phi();
      nEvent.genEta1    = sortedGLeptons[i1].p.Eta();
      nEvent.genMID1gen = fTR->GenLeptonMID[sortedGLeptons[i1].index];
      nEvent.genGMID1gen = fTR->GenLeptonGMID[sortedGLeptons[i1].index];
      if(i2<sortedGLeptons.size())
        {
          TLorentzVector genZvector = sortedGLeptons[i1].p + sortedGLeptons[i2].p;
          nEvent.genRecoil  = (-GenMETvector - genZvector).Pt();
          nEvent.genPt2     = sortedGLeptons[i2].p.Pt();
          nEvent.genId2     = sortedGLeptons[i2].type;
          nEvent.genPhi2    = sortedGLeptons[i2].p.Phi();
          nEvent.genEta2    = sortedGLeptons[i2].p.Eta();
          nEvent.genMID2gen = fTR->GenLeptonMID[sortedGLeptons[i2].index];
          nEvent.genGMID2gen= fTR->GenLeptonGMID[sortedGLeptons[i1].index];
          nEvent.genZPt     = genZvector.Pt();
          nEvent.genDRll    = sortedGLeptons[i2].p.DeltaR(sortedGLeptons[i1].p);
          nEvent.genMll     = genZvector.M();
          nEvent.genJZB     = nEvent.genRecoil - genZvector.Pt();
	  nEvent.dphigenZgenMet = (sortedGLeptons[i1].p + sortedGLeptons[i2].p).DeltaPhi(GenMETvector);
        }

      if(sortedGLeptons.size()>2) {
	int genTriLepton1 = 0;
	int genTriLepton2 = 1;
	int genTriLepton3 = 2;
	
	int genBadTriLepton1 = 0;
	int genBadTriLepton2 = 1;
	int genBadTriLepton3 = 1;
	
	float mindiff=-1;
	int BestCandidate=sortedGLeptons.size();
	// Trileptons
	for(; genTriLepton2 < sortedGLeptons.size(); genTriLepton2++) {
	  if(sortedGLeptons[0].charge*sortedGLeptons[genTriLepton2].charge<0 &&  sortedGLeptons[0].type==sortedGLeptons[genTriLepton2].type ) {
	    float curr_mindiff=abs((sortedGLeptons[0].p+sortedGLeptons[genTriLepton2].p).M()-91.2);
	    if(mindiff<0||curr_mindiff<mindiff) {
	      mindiff=curr_mindiff;
	      BestCandidate=genTriLepton2;
	    }
	  }
	}
	
	genTriLepton2=BestCandidate;
	
	
	for(; genTriLepton3 < sortedGLeptons.size(); genTriLepton3++) {
	  if(genTriLepton3==genTriLepton1 || genTriLepton3==genTriLepton2) continue;
	  else break; // pick the first lepton that's left!
	}
  
	mindiff=-1;
	BestCandidate=sortedGLeptons.size();
	// Trileptons
	for(; genBadTriLepton2 < sortedGLeptons.size(); genBadTriLepton2++) {
	  if(sortedGLeptons[0].charge*sortedGLeptons[genBadTriLepton2].charge<0) {
	    float curr_mindiff=abs((sortedGLeptons[0].p+sortedGLeptons[genBadTriLepton2].p).M()-91.2);
	    if(mindiff<0||curr_mindiff<mindiff) {
	      mindiff=curr_mindiff;
	      BestCandidate=genBadTriLepton2;
	    }
	  }
	}
	
	genBadTriLepton2=BestCandidate;
	
	
	for(; genBadTriLepton3 < sortedGLeptons.size(); genBadTriLepton3++) {
	  if(genBadTriLepton3==genBadTriLepton1 || genBadTriLepton3==genBadTriLepton2) continue;
	  else break; // pick the first lepton that's left!
	}
	
	if(sortedGLeptons.size()>=3 && genTriLepton1<sortedGLeptons.size() && genTriLepton2<sortedGLeptons.size() && genTriLepton3<sortedGLeptons.size() ) {
	  TLorentzVector lp1(sortedGLeptons[genTriLepton1].p);
	  TLorentzVector lp2(sortedGLeptons[genTriLepton2].p);
	  TLorentzVector lp3(sortedGLeptons[genTriLepton3].p);
	  
	  nEvent.gentri_pt1=sortedGLeptons[genTriLepton1].p.Pt();
	  nEvent.gentri_pt2=sortedGLeptons[genTriLepton2].p.Pt();
	  nEvent.gentri_pt3=sortedGLeptons[genTriLepton3].p.Pt();
	  
	  nEvent.gentri_id1=sortedGLeptons[genTriLepton1].type;
	  nEvent.gentri_id2=sortedGLeptons[genTriLepton2].type;
	  nEvent.gentri_id3=sortedGLeptons[genTriLepton3].type;
	  
	  nEvent.gentri_ch1=sortedGLeptons[genTriLepton1].charge;
	  nEvent.gentri_ch2=sortedGLeptons[genTriLepton2].charge;
	  nEvent.gentri_ch3=sortedGLeptons[genTriLepton3].charge;
	  
	  nEvent.gentri_mlll=(lp1+lp2+lp3).M();
	  nEvent.gentri_mll=(lp1+lp2).M();
	  nEvent.gentri_submll=(lp2+lp3).M();
	  
	  nEvent.gentri_badmll=(sortedGLeptons[genBadTriLepton1].p+sortedGLeptons[genBadTriLepton2].p).M();
	  nEvent.gentri_badsubmll=(sortedGLeptons[genBadTriLepton2].p+sortedGLeptons[genBadTriLepton3].p).M();
    
	  float AngleBetweenMETandThirdLepton=sortedGLeptons[genTriLepton3].p.DeltaPhi(GenMETvector);
	  nEvent.gentri_mT = sqrt(2 * nEvent.genMET * sortedGLeptons[genTriLepton3].p.Pt() * ( 1 - cos(AngleBetweenMETandThirdLepton)));
	  
	  AngleBetweenMETandThirdLepton=sortedGLeptons[genBadTriLepton3].p.DeltaPhi(GenMETvector);
	  nEvent.gentri_badmT = sqrt(2 * nEvent.genMET * sortedGLeptons[genBadTriLepton3].p.Pt() * ( 1 - cos(AngleBetweenMETandThirdLepton)));
	  
	  nEvent.gentri_genMID1 = fTR->GenLeptonMID[sortedGLeptons[genTriLepton1].index];
	  nEvent.gentri_genMID2 = fTR->GenLeptonMID[sortedGLeptons[genTriLepton2].index];
	  nEvent.gentri_genMID3 = fTR->GenLeptonMID[sortedGLeptons[genTriLepton3].index];
	  
	  nEvent.gentri_badgenMID1 = fTR->GenLeptonMID[sortedGLeptons[genBadTriLepton1].index];
	  nEvent.gentri_badgenMID2 = fTR->GenLeptonMID[sortedGLeptons[genBadTriLepton2].index];
	  nEvent.gentri_badgenMID3 = fTR->GenLeptonMID[sortedGLeptons[genBadTriLepton3].index];
	  
	  if(abs(nEvent.gentri_genMID1)==15) nEvent.gentri_genMID1=fTR->GenLeptonGMID[sortedGLeptons[genTriLepton1].index];
	  if(abs(nEvent.gentri_genMID2)==15) nEvent.gentri_genMID2=fTR->GenLeptonGMID[sortedGLeptons[genTriLepton2].index];
	  if(abs(nEvent.gentri_genMID3)==15) nEvent.gentri_genMID3=fTR->GenLeptonGMID[sortedGLeptons[genTriLepton3].index];
	  
	  if(abs(nEvent.gentri_badgenMID1)==15) nEvent.gentri_badgenMID1=fTR->GenLeptonGMID[sortedGLeptons[genBadTriLepton1].index];
	  if(abs(nEvent.gentri_badgenMID2)==15) nEvent.gentri_badgenMID2=fTR->GenLeptonGMID[sortedGLeptons[genBadTriLepton2].index];
	  if(abs(nEvent.gentri_badgenMID3)==15) nEvent.gentri_badgenMID3=fTR->GenLeptonGMID[sortedGLeptons[genBadTriLepton3].index];
	  
	  nEvent.gentri_GoodZMatch=false;
	  nEvent.gentri_GoodWMatch=false;
	  if(nEvent.gentri_genMID1==23 && nEvent.gentri_genMID1 ==23) nEvent.gentri_GoodZMatch=true;
	  if(abs(nEvent.gentri_genMID3)==24 ) nEvent.gentri_GoodWMatch=true;
	    
	  nEvent.gentri_dR12 = lp1.DeltaR(lp2);
	  nEvent.gentri_dR13 = lp1.DeltaR(lp3); 
	  nEvent.gentri_dR23 = lp2.DeltaR(lp3);
	  
	  nEvent.gentri_eta1 = lp1.Eta();
	  nEvent.gentri_eta2 = lp2.Eta();
	  nEvent.gentri_eta3 = lp3.Eta();
	  
	} // end of three or more usable leptons
	  // ** and now the pure generator one (i.e. no combinatorial errors!) -- only for WZ otherwise this is pointless!
	  
	int pgenTriLepton1 = 0;
	int pgenTriLepton2 = 0;
	int pgenTriLepton3 = 0;
	
	for(; pgenTriLepton1 < sortedGLeptons.size(); pgenTriLepton1++) {
	  int mother=abs(fTR->GenLeptonMID[sortedGLeptons[pgenTriLepton1].index]);
	  if(mother==15) mother=abs(fTR->GenLeptonGMID[sortedGLeptons[pgenTriLepton1].index]);
	  if(mother==23) break; // found a lepton coming from a Z
	}
	for(pgenTriLepton2=pgenTriLepton1+1; pgenTriLepton2 < sortedGLeptons.size(); pgenTriLepton2++) {
	  int mother=abs(fTR->GenLeptonMID[sortedGLeptons[pgenTriLepton2].index]);
	  if(mother==15) mother=abs(fTR->GenLeptonGMID[sortedGLeptons[pgenTriLepton2].index]);
	  if(mother==23) break; // found a lepton coming from a Z
	}
	for(; pgenTriLepton3 < sortedGLeptons.size(); pgenTriLepton3++) {
	  int mother=abs(fTR->GenLeptonMID[sortedGLeptons[pgenTriLepton3].index]);
	  if(mother==15) mother=abs(fTR->GenLeptonGMID[sortedGLeptons[pgenTriLepton3].index]);
	  if(mother==24) break; // found a lepton coming from a W
	}
	
	if(sortedGLeptons.size()>=3 && pgenTriLepton1<sortedGLeptons.size() && pgenTriLepton2<sortedGLeptons.size() && pgenTriLepton3<sortedGLeptons.size() ) {
	  TLorentzVector lp1(sortedGLeptons[pgenTriLepton1].p);
	  TLorentzVector lp2(sortedGLeptons[pgenTriLepton2].p);
	  TLorentzVector lp3(sortedGLeptons[pgenTriLepton3].p);
	  
	  nEvent.pgentri_pt1=sortedGLeptons[pgenTriLepton1].p.Pt();
	  nEvent.pgentri_pt2=sortedGLeptons[pgenTriLepton2].p.Pt();
	  nEvent.pgentri_pt3=sortedGLeptons[pgenTriLepton3].p.Pt();
	  
	  nEvent.pgentri_id1=sortedGLeptons[pgenTriLepton1].type;
	  nEvent.pgentri_id2=sortedGLeptons[pgenTriLepton2].type;
	  nEvent.pgentri_id3=sortedGLeptons[pgenTriLepton3].type;
	  
	  nEvent.pgentri_ch1=sortedGLeptons[pgenTriLepton1].charge;
	  nEvent.pgentri_ch2=sortedGLeptons[pgenTriLepton2].charge;
	  nEvent.pgentri_ch3=sortedGLeptons[pgenTriLepton3].charge;
	  
	  nEvent.pgentri_mlll=(lp1+lp2+lp3).M();
	  nEvent.pgentri_mll=(lp1+lp2).M();
	  nEvent.pgentri_submll=(lp2+lp3).M();
	  
	  float AngleBetweenMETandThirdLepton=sortedGLeptons[pgenTriLepton3].p.DeltaPhi(GenMETvector);
	  nEvent.pgentri_mT = sqrt(2 * nEvent.genMET * sortedGLeptons[pgenTriLepton3].p.Pt() * ( 1 - cos(AngleBetweenMETandThirdLepton)));
	  
	  nEvent.pgentri_dR12 = lp1.DeltaR(lp2);
	  nEvent.pgentri_dR13 = lp1.DeltaR(lp3); 
	  nEvent.pgentri_dR23 = lp2.DeltaR(lp3);
	  
	  nEvent.pgentri_eta1 = lp1.Eta();
	  nEvent.pgentri_eta2 = lp2.Eta();
	  nEvent.pgentri_eta3 = lp3.Eta();
	  
	} // end of Z/W pure gen info
      }//end of if there are three or more gleptons
    }//end of if there are any gleptons
}
     
    
int JZBAnalysis::DetermineFlavor(bool fdoGenInfo,TreeReader *fTR) {
  int flavorCounter[7];
  for(int i=0;i<7;i++) flavorCounter[i]=0;
  for(int i=0;i<fTR->nGenParticles;i++) {
    if(fTR->genInfoStatus[i]!=3) continue;
    if(abs(fTR->genInfoId[i])>10||abs(fTR->genInfoId[i])==0) continue; // not a quark
    flavorCounter[abs(fTR->genInfoId[i])]++;
  }
  for(int i=6;i>0;i--) {
    if(flavorCounter[i]>0) return i;
  }
  return 1;
}

bool JZBAnalysis::DecaysToTaus(bool fdoGenInfo,TreeReader *fTR) {
  //built in potential to look at lepton flavors other than taus
  int nTaus=0;
  int nNeutrinos=0;
  int nHad=0;
  int nMu=0;
  int nEl=0;
  for(int iMother=0;iMother<fTR->nGenParticles;iMother++) {
    if(fTR->genInfoStatus[iMother]==2) continue;
    if(!(abs(fTR->genInfoId[iMother])==23)) continue; // not a Z
    for(int da=iMother+1;da<fTR->nGenParticles;da++) {
      if(fTR->genInfoStatus[da]==2) continue;
      if(!(fTR->genInfoMo1[da]==iMother) ) continue; // not a daughter
      int Pid=abs(fTR->genInfoId[da]);
      if(Pid==11) nEl++;
      if(Pid==13) nMu++;
      if(Pid==15) nTaus++;
      if(Pid==12) nNeutrinos++;//e neutrino
      if(Pid==14) nNeutrinos++;//m neutrino
      if(Pid==16) nNeutrinos++;//t neutrino
      if(Pid>0&&Pid<10) nHad++;
    }//end of daughter loop
  }//end of particle loop
  if(nTaus>0) return true;
  return false;
}
  


bool JZBAnalysis::MatchTrigger(lepton *a) {

  //cout << "We are testing lepton: " << a->p << " which is of type: " << a->type << endl;
  if(a->type == 0) {
    if(! passTriggers(singleElTriggerPaths)) return false; 
    for(int k = 0; k < fTR->NHLTObjs[1]; k++) {
      if(abs(fTR->HLTObjectID1[k]) == 11) {
        TLorentzVector tmpVector(0, 0, 0, 0);
        tmpVector.SetPtEtaPhiM(fTR->HLTObjectPt1[k], fTR->HLTObjectEta1[k],
                               fTR->HLTObjectPhi1[k], a->p.M());
        if(tmpVector.DeltaR(a->p)<0.1) {
          return true;
        }
      }
    }
  } else {
    if(! passTriggers(singleMuTriggerPaths)) return false; 
    for(int k = 0; k < fTR->NHLTObjs[0]; k++) {
      if(abs(fTR->HLTObjectID0[k]) == 13) {
        TLorentzVector tmpVector(0, 0, 0, 0);
        tmpVector.SetPtEtaPhiM(fTR->HLTObjectPt0[k], fTR->HLTObjectEta0[k],
                               fTR->HLTObjectPhi0[k], a->p.M());
        if(tmpVector.DeltaR(a->p)<0.1) {
          return true;
        }
      }
    }
  }
  return false;

}



