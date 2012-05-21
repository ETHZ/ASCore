#include "helper/Utilities.hh"
#include "JZBAnalysis.hh"
#include "TF1.h"
#include <time.h>
#include <TRandom.h>
#include "TF1.h"
#include "TTree.h"
#include <cstdlib>
#include <assert.h>

using namespace std;


#define jMax 30  // do not touch this
#define rMax 30
#define Zmax 30

enum METTYPE { mettype_min, RAW = mettype_min, DUM, TCMET, MUJESCORRMET, PFMET, SUMET, PFRECOILMET, RECOILMET, mettype_max };
enum JZBTYPE { jzbtype_min, CALOJZB = jzbtype_min, PFJZB, RECOILJZB, PFRECOILJZB, TCJZB, jzbtype_max };

string sjzbversion="$Revision: 1.70.2.20 $";
string sjzbinfo="";

float firstLeptonPtCut  = 10.0;
float secondLeptonPtCut = 10.0;

/*
$Id$
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
  float pt;
  float phi;
  bool is_data;

  float pt1; // leading leptons
  float pt2;
  float iso1;
  float iso2;
  bool  isConv1; // Photon conversion flag
  bool  isConv2;

  float genPt1; // leading legenPtons
  float genPt2;
  int   genId1;
  int   genId2;
  int   genMID;
  int   genMID1;
  int   genMID2;
  int   genMID3;
  int   genMID4;
  int   genMID5;
  int   genGMID;
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
  int   genNjetsTwoSix;
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
  float dphipfRecoilMet1;
  float dphipfRecoilMet2;

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
  float pfHT;
  float pfGoodHT;
  float pfTightHT;

  int pfJetGoodNum;
  int pfJetGoodNumID;
  int pfJetGoodNump1sigma;
  int pfJetGoodNumn1sigma;
  float pfJetGoodPt[jMax];
  float pfJetGoodEta[jMax];
  float pfJetGoodPhi[jMax];
  bool  pfJetGoodID[jMax];
  float bTagProbCSVBP[jMax];
  float bTagProbCSVMVA[jMax];
  
  int pfJetGoodNumBtag;
  int pfJetGoodNumIDBtag;
  float pfJetGoodPtBtag[jMax];
  float pfJetGoodEtaBtag[jMax];
  float pfJetGoodPhiBtag[jMax];
  bool  pfJetGoodIDBtag[jMax];

  int pfJetGoodNum40;
  int pfJetGoodNum50;
  int pfJetGoodNum60;


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
  float PUweight;
  bool passed_triggers;
  int trigger_bit;
  float mGlu;
  float mChi;
  float mLSP;
  float xSMS;
  float xbarSMS;
  float mGMSBGlu;
  float mGMSBChi;
  float mGMSBLSP;
  float M0;
  float A0;
  float M12;
  float signMu;
  
  
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

};

nanoEvent::nanoEvent(){};
void nanoEvent::reset()
{

  mll=0; // di-lepton system
  pt=0;
  phi=0;

  is_data=false;
  NPdfs=0;
  pdfWsum=0;

  process=0;

  pt1=0;
  pt2=0;
  iso1=0;
  iso2=0;
  isConv1 = false;
  isConv2 = false;
  genPt1=0;
  genPt2=0;
  genEta1=0;
  genEta2=0;
  genId1=0;
  genId2=0;
  genMID=0;
  genMID1=0;
  genMID2=0;
  genMID3=0;
  genMID4=0;
  genMID5=0;
  genGMID=0;
  genGMID1=0;
  genGMID2=0;
  genMET=0;
  genZPt=0;
  genMll=0;
  genRecoil=0;
  genJZB = 0;
  genNjets = 0;
  genNjetsTwoSix = 0;
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
  }
  pfJetNum=0;
  pfHT=0;
  pfGoodHT=0;
  pfTightHT=0;

  for(int jCounter=0;jCounter<jMax;jCounter++){
    pfJetGoodPt[jCounter]=0;
    pfJetGoodEta[jCounter]=0;
    pfJetGoodPhi[jCounter]=0;
    pfJetGoodID[jCounter]=0;
    pfJetGoodPtBtag[jCounter]=0;
    pfJetGoodEtaBtag[jCounter]=0;
    pfJetGoodPhiBtag[jCounter]=0;
    pfJetGoodIDBtag[jCounter]=0;
    bTagProbCSVBP[jCounter] = 0;
    bTagProbCSVMVA[jCounter] = 0;
  }

  pfJetGoodNum=0;
  pfJetGoodNumID=0;
  pfJetGoodNump1sigma=0;
  pfJetGoodNumn1sigma=0;

  pfJetGoodNumBtag=0;
  pfJetGoodNum40=0;
  pfJetGoodNum50=0;
  pfJetGoodNum60=0;

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
  Efficiencyweightonly = 1.0;
  weightEffDown = 1.0;
  weightEffUp = 1.0;

  mGlu=0;
  mChi=0;
  mLSP=0;
  xSMS=0;
  xbarSMS=0;
  mGMSBGlu=0;
  mGMSBChi=0;
  mGMSBLSP=0;
  A0=0;
  M0=0;
  M12=0;
  signMu=0;
  
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

}


TTree *myTree;
TTree *FullTree;
TTree *myInfo;

nanoEvent nEvent;


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


JZBAnalysis::JZBAnalysis(TreeReader *tr, std::string dataType, bool fullCleaning, bool isModelScan, bool makeSmall, bool doGenInfo) :
  UserAnalysisBase(tr), fDataType_(dataType), fFullCleaning_(fullCleaning) , fisModelScan(isModelScan) , fmakeSmall(makeSmall), fdoGenInfo(doGenInfo) {
  // Define trigger paths to check
  addPath(elTriggerPaths,"HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL",1,8);
  addPath(elTriggerPaths,"HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL",1,5);
  addPath(elTriggerPaths,"HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL",1,10);
  addPath(elTriggerPaths,"HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL",11, 20);

  addPath(muTriggerPaths,"HLT_DoubleMu6",1,8); // v5 is actually not used
  addPath(muTriggerPaths,"HLT_DoubleMu7",1,2); // v1,2,5,7,8,11,12 are really used
  addPath(muTriggerPaths,"HLT_DoubleMu7",5,8); 
  addPath(muTriggerPaths,"HLT_DoubleMu7",11,12);
  addPath(muTriggerPaths,"HLT_Mu13_Mu8",2,4); // 2,3,4,6,7,10,11
  addPath(muTriggerPaths,"HLT_Mu13_Mu8",6,7);
  addPath(muTriggerPaths,"HLT_Mu13_Mu8",10,11);
  addPath(muTriggerPaths,"HLT_Mu17_Mu8",2,4); // 2,3,4,6,7,10,11
  addPath(muTriggerPaths,"HLT_Mu17_Mu8",6,7);
  addPath(muTriggerPaths,"HLT_Mu17_Mu8",10,11);
  addPath(muTriggerPaths,"HLT_Mu17_Mu8",12,16);
  
  addPath(emTriggerPaths,"HLT_Mu17_Ele8_CaloIdL",1,9);
  addPath(emTriggerPaths,"HLT_Mu17_Ele8_CaloIdL",12,13);
  addPath(emTriggerPaths,"HLT_Mu8_Ele17_CaloIdL",1,9);
  addPath(emTriggerPaths,"HLT_Mu8_Ele17_CaloIdL",12,13);
  addPath(emTriggerPaths,"HLT_Mu8_Ele17_CaloIdT_CaloIsoVL",1,1);
  addPath(emTriggerPaths,"HLT_Mu8_Ele17_CaloIdT_CaloIsoVL",3,3);
  addPath(emTriggerPaths,"HLT_Mu8_Ele17_CaloIdT_CaloIsoVL",4,4);
  addPath(emTriggerPaths,"HLT_Mu8_Ele17_CaloIdT_CaloIsoVL",7,8);
  addPath(emTriggerPaths,"HLT_Mu17_Ele8_CaloIdT_CaloIsoVL",1,1);
  addPath(emTriggerPaths,"HLT_Mu17_Ele8_CaloIdT_CaloIsoVL",3,3);
  addPath(emTriggerPaths,"HLT_Mu17_Ele8_CaloIdT_CaloIsoVL",4,4);
  addPath(emTriggerPaths,"HLT_Mu17_Ele8_CaloIdT_CaloIsoVL",7,8);
  addPath(emTriggerPaths,"HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL", 1, 10);
  addPath(emTriggerPaths,"HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL", 1, 10);
  
}

//________________________________________________________________________________________
JZBAnalysis::~JZBAnalysis(){}

//________________________________________________________________________________________
void JZBAnalysis::Begin(TFile *f){

  f->cd();
  rand_ = new TRandom();

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
  
  myInfo->Fill();
  myInfo->Write();

  FullTree = new TTree("Allevents","Allevents");
  FullTree->Branch("is_data",&nEvent.is_data,"is_data/O");
  
  myTree = new TTree("events","events");

  myTree->Branch("is_data",&nEvent.is_data,"is_data/O");
  myTree->Branch("mll",&nEvent.mll,"mll/F");
  myTree->Branch("pt",&nEvent.pt,"pt/F");
  myTree->Branch("phi",&nEvent.phi,"phi/F");
  myTree->Branch("pt1",&nEvent.pt1,"pt1/F");
  myTree->Branch("pt2",&nEvent.pt2,"pt2/F");
  myTree->Branch("iso1",&nEvent.iso1,"iso1/F");
  myTree->Branch("iso2",&nEvent.iso2,"iso2/F");
  myTree->Branch("isConv1",&nEvent.isConv1,"isConv1/O");
  myTree->Branch("isConv2",&nEvent.isConv2,"isConv2/O");
  myTree->Branch("genPt1",&nEvent.genPt1,"genPt1/F");
  myTree->Branch("genPt2",&nEvent.genPt2,"genPt2/F");
  myTree->Branch("genEta1",&nEvent.genEta1,"genEta1/F");
  myTree->Branch("genEta2",&nEvent.genEta2,"genEta2/F");
  myTree->Branch("genId1",&nEvent.genId1,"genId1/I");
  myTree->Branch("genId2",&nEvent.genId2,"genId2/I");
  myTree->Branch("genMID",&nEvent.genMID,"genMID/I");
  myTree->Branch("genMID1",&nEvent.genMID1,"genMID1/I");
  myTree->Branch("genMID2",&nEvent.genMID2,"genMID2/I");
  myTree->Branch("genMID3",&nEvent.genMID3,"genMID3/I");
  myTree->Branch("genMID4",&nEvent.genMID4,"genMID4/I");
  myTree->Branch("genMID5",&nEvent.genMID5,"genMID5/I");
  myTree->Branch("genGMID",&nEvent.genGMID,"genGMID/I");
  myTree->Branch("genGMID1",&nEvent.genGMID1,"genGMID1/I");
  myTree->Branch("genGMID2",&nEvent.genGMID2,"genGMID2/I");
  myTree->Branch("genMET",&nEvent.genMET,"genMET/F");
  myTree->Branch("genZPt",&nEvent.genZPt,"genZPt/F");
  myTree->Branch("genMll",&nEvent.genMll,"genMll/F");
  myTree->Branch("genRecoil",&nEvent.genRecoil,"genRecoil/F");
  myTree->Branch("genJZB",&nEvent.genJZB,"genJZB/F");
  myTree->Branch("genNjets",&nEvent.genNjets,"genNjets/I");
  myTree->Branch("genNjetsTwoSix",&nEvent.genNjetsTwoSix,"genNjetsTwoSix/I");
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
  myTree->Branch("pfHT",&nEvent.pfHT,"pfHT/F");
  myTree->Branch("pfGoodHT",&nEvent.pfGoodHT,"pfGoodHT/F");
  myTree->Branch("pfTightHT",&nEvent.pfTightHT,"pfTightHT/F");

  myTree->Branch("pfJetGoodNum",&nEvent.pfJetGoodNum,"pfJetGoodNum/I");
  myTree->Branch("pfJetGoodNumBtag",&nEvent.pfJetGoodNumBtag,"pfJetGoodNumBtag/I");
  myTree->Branch("pfJetGoodNumID",&nEvent.pfJetGoodNumID,"pfJetGoodNumID/I");
  myTree->Branch("pfJetGoodNump1sigma",&nEvent.pfJetGoodNump1sigma,"pfJetGoodNump1sigma/I");
  myTree->Branch("pfJetGoodNumn1sigma",&nEvent.pfJetGoodNumn1sigma,"pfJetGoodNumn1sigma/I");

  myTree->Branch("pfJetGoodPt", nEvent.pfJetGoodPt,"pfJetGoodPt[pfJetGoodNum]/F");
  myTree->Branch("pfJetGoodEta",nEvent.pfJetGoodEta,"pfJetGoodEta[pfJetGoodNum]/F");
  myTree->Branch("pfJetGoodPhi",nEvent.pfJetGoodPhi,"pfJetGoodPhi[pfJetGoodNum]/F");
  myTree->Branch("pfJetGoodID", nEvent.pfJetGoodID,"pfJetGoodID[pfJetGoodNum]/O");

  myTree->Branch("pfJetGoodPtBtag", nEvent.pfJetGoodPtBtag,"pfJetGoodPtBag[pfJetGoodNumBtag]/F");
  myTree->Branch("pfJetGoodEtaBtag",nEvent.pfJetGoodEtaBtag,"pfJetGoodEtaBtag[pfJetGoodNumBtag]/F");
  myTree->Branch("pfJetGoodPhiBtag",nEvent.pfJetGoodPhiBtag,"pfJetGoodPhiBtag[pfJetGoodNumBtag]/F");
  myTree->Branch("pfJetGoodIDBtag", nEvent.pfJetGoodIDBtag,"pfJetGoodID[pfJetGoodNumBtag]/O");
  myTree->Branch("bTagProbCSVBP", nEvent.bTagProbCSVBP,"bTagProbCSVBP[pfJetGoodNum]/F");
  myTree->Branch("bTagProbCSVMVA", nEvent.bTagProbCSVMVA,"bTagProbCSVMVA[pfJetGoodNum]/F");

  myTree->Branch("pfJetGoodNum40",&nEvent.pfJetGoodNum40,"pfJetGoodNum40/I");
  myTree->Branch("pfJetGoodNum50",&nEvent.pfJetGoodNum50,"pfJetGoodNum50/I");
  myTree->Branch("pfJetGoodNum60",&nEvent.pfJetGoodNum60,"pfJetGoodNum60/I");

  myTree->Branch("jzb",nEvent.jzb,Form("jzb[%d]/F",int(jzbtype_max)));
  myTree->Branch("sjzb",nEvent.sjzb,Form("sjzb[%d]/F",int(jzbtype_max)));
  myTree->Branch("dphi_sumJetVSZ",nEvent.dphi_sumJetVSZ,Form("dphi_sumJetVSZ[%d]/F",int(jzbtype_max)));
  myTree->Branch("sumJetPt",nEvent.sumJetPt,Form("sumJetPt[%d]/F",int(jzbtype_max)));

  myTree->Branch("weight", &nEvent.weight,"weight/F");
  myTree->Branch("PUweight",&nEvent.PUweight,"PUweight/F");
  myTree->Branch("Efficiencyweightonly",&nEvent.Efficiencyweightonly,"Efficiencyweightonly/F");
  myTree->Branch("weightEffDown",&nEvent.weightEffDown,"weightEffDown/F");
  myTree->Branch("weightEffUp",&nEvent.weightEffUp,"weightEffUp/F");

  myTree->Branch("passed_triggers", &nEvent.passed_triggers,"passed_triggers/O");
  myTree->Branch("trigger_bit", &nEvent.trigger_bit,"trigger_bit/I");
  myTree->Branch("MassGlu",&nEvent.mGlu,"MassGlu/F");
  myTree->Branch("MassChi",&nEvent.mChi,"MassChi/F");
  myTree->Branch("MassLSP",&nEvent.mLSP,"MassLSP/F");
  myTree->Branch("xSMS",&nEvent.xSMS,"xSMS/F");
  myTree->Branch("xbarSMS",&nEvent.xbarSMS,"xbarSMS/F");
  myTree->Branch("MassGMSBGlu",&nEvent.mGMSBGlu,"MassGlu/F");
  myTree->Branch("MassGMSBChi",&nEvent.mGMSBChi,"MassChi/F");
  myTree->Branch("MassGMSBLSP",&nEvent.mGMSBLSP,"MassLSP/F");
  myTree->Branch("M0",&nEvent.M0,"M0/F");
  myTree->Branch("A0",&nEvent.A0,"A0/F");
  myTree->Branch("M12",&nEvent.M12,"M12/F");
  myTree->Branch("signMu",&nEvent.signMu,"signMu/F");
  myTree->Branch("NPdfs",&nEvent.NPdfs,"NPdfs/I");
  myTree->Branch("pdfW",nEvent.pdfW,"pdfW[NPdfs]/F");
  myTree->Branch("pdfWsum",&nEvent.pdfWsum,"pdfWsum/F");
  
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

  counters[EV].setName("Events");
  counters[TR].setName("Triggers");
  counters[MU].setName("Muons");
  counters[EL].setName("Electrons");
  counters[PJ].setName("PFJets");

  // Define counters (so we have them in the right order)
  counters[EV].fill("All events",0.);
  if ( fDataType_ != "mc" ) {
    counters[EV].fill("... pass electron triggers",0.);
    counters[EV].fill("... pass muon triggers",0.);
    counters[EV].fill("... pass EM triggers",0.);
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
  assert(foundUnprescaled);

  return passed;

}

bool is_neutrino(int code) {
  if(abs(code)==12) return true; // electron neutrino
  if(abs(code)==14) return true; // muon neutrino
  if(abs(code)==16) return true; // tau neutrino
  if(abs(code)==18) return true; // tau' neutrino
  return false;
}

bool is_charged_lepton(int code) {
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
  FullTree->Fill();


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
        nEvent.A0=fTR->A0;
        nEvent.M0=fTR->M0;
        nEvent.signMu=fTR->signMu;
        nEvent.M12=fTR->M12;
        nEvent.NPdfs=fTR->NPdfs;
        for(int i=0;i<fTR->NPdfs;i++) nEvent.pdfW[i]=fTR->pdfW[i];
	nEvent.pdfWsum=fTR->pdfWsum; 
      } else {
	//don't attempt to do PURW for model scans
	nEvent.PUweight  = GetPUWeight(fTR->PUnumInteractions);
	nEvent.weight    = GetPUWeight(fTR->PUnumInteractions);
      }
      
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
	for(int i=0;i<nGenParticles&&fdoGenInfo;i++) {
	  if(fTR->genInfoStatus[i]!=3) continue;
	  int thisParticleId = fTR->genInfoId[i];
	  if(fdoGenInfo&&abs(thisParticleId)==23) {
	    //dealing with a Z
	    int motherIndex=fTR->genInfoMo1[i];
	    if(motherIndex>=0) nEvent.SourceOfZ[nEvent.nZ]=fTR->genInfoId[motherIndex];
	    nEvent.nZ++;
	    for(int da=0;da<fTR->nGenParticles;da++) {
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

	if(fdoGenInfo) {
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
	}//end of if(f_doGenInfo)

	
	if(nchimass>0&&nlspmass>0&&nglumass>0)  nEvent.realx=(chimass/nchimass - lspmass/nlspmass)/(glumass/nglumass-lspmass/nlspmass);
	//at this point we use the fact that one of the three bits of information in the LHE event comment is the imposed x - the only question is which one ;-) 
	//note: the bit of information in the comment is actually xbar, so we need to store 1-xbar to get our definition of x. 
	if(nEvent.mGlu>0 && nEvent.mGlu<1) nEvent.imposedx=1-nEvent.mGlu;
	if(nEvent.mChi>0 && nEvent.mChi<1) nEvent.imposedx=1-nEvent.mChi;
	if(nEvent.mLSP>0 && nEvent.mLSP<1) nEvent.imposedx=1-nEvent.mLSP;
  } // end of mc if
  
  // Trigger information
  nEvent.passed_triggers=0;
  if ( fDataType_ != "mc" ) 
    {
      nEvent.is_data=true;
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
    }

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
    
  // Check for OS combination
  for(; PosLepton2 < sortedGoodLeptons.size(); PosLepton2++) {
    if(sortedGoodLeptons[0].charge*sortedGoodLeptons[PosLepton2].charge<0) break;
  }
  if(PosLepton2 == sortedGoodLeptons.size()) {
    if (isMC&&!fmakeSmall) myTree->Fill();
    return;
  }
  counters[EV].fill("... has at least 2 OS leptons");

  // Preselection
  if(sortedGoodLeptons[PosLepton1].p.Pt() > firstLeptonPtCut && sortedGoodLeptons[PosLepton2].p.Pt() > secondLeptonPtCut) {

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
    
    nEvent.mll=(sortedGoodLeptons[PosLepton2].p+sortedGoodLeptons[PosLepton1].p).M();
    nEvent.phi=(sortedGoodLeptons[PosLepton2].p+sortedGoodLeptons[PosLepton1].p).Phi();
    nEvent.pt=(sortedGoodLeptons[PosLepton2].p+sortedGoodLeptons[PosLepton1].p).Pt();
    nEvent.dphi=sortedGoodLeptons[PosLepton2].p.DeltaPhi(sortedGoodLeptons[PosLepton1].p);
    
    nEvent.ElCInfoIsGsfCtfCons=sortedGoodLeptons[PosLepton2].ElCInfoIsGsfCtfCons&&sortedGoodLeptons[PosLepton1].ElCInfoIsGsfCtfCons;
    nEvent.ElCInfoIsGsfCtfScPixCons=sortedGoodLeptons[PosLepton2].ElCInfoIsGsfCtfScPixCons&&sortedGoodLeptons[PosLepton1].ElCInfoIsGsfCtfScPixCons;
    nEvent.ElCInfoIsGsfScPixCons=sortedGoodLeptons[PosLepton2].ElCInfoIsGsfScPixCons&&sortedGoodLeptons[PosLepton1].ElCInfoIsGsfScPixCons;

    float lepweightErr;
    float lepweight=GetLeptonWeight(nEvent.id1,nEvent.phi1,nEvent.eta1,nEvent.id2,nEvent.phi2,nEvent.eta2,lepweightErr);
    
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
  
  float tcMETpx = fTR->PFType1METpx;
  float tcMETpy = fTR->PFType1METpy;
  
  TLorentzVector pfMETvector(pfMETpx,pfMETpy,0,0);
  TLorentzVector type1METvector(tcMETpx,tcMETpy,0,0);
  TLorentzVector sumOfPFJets(0,0,0,0);
  nEvent.pfJetNum=0;
  nEvent.pfJetGoodNum=0;
  nEvent.pfJetGoodNumBtag=0;
  nEvent.pfJetGoodNum40=0;
  nEvent.pfJetGoodNum50=0;
  nEvent.pfJetGoodNum60=0;
  
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
      bool  isJetID = IsGoodBasicPFJet(i,false,3.0);
      
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
      
      // Keep jets over min. pt threshold
      if ( !(jpt>20) ) continue;
      counters[PJ].fill("... pt>20.");

      //Get Uncertainty
      jecUnc->setJetEta(jeta);
      jecUnc->setJetPt(jpt); // here you must use the CORRECTED jet pt
      float unc = jecUnc->getUncertainty(true); 
      nEvent.pfJetPt[nEvent.pfJetNum]    = jpt;
      nEvent.pfJetEta[nEvent.pfJetNum]   = jeta;
      nEvent.pfJetPhi[nEvent.pfJetNum]   = jphi;
      nEvent.pfJetScale[nEvent.pfJetNum] = jesC;
      nEvent.pfJetScaleUnc[nEvent.pfJetNum] = unc;
      nEvent.pfJetID[nEvent.pfJetNum]    = isJetID;
      nEvent.pfJetDphiMet[nEvent.pfJetNum] = aJet.DeltaPhi(pfMETvector);
      nEvent.pfJetNum = nEvent.pfJetNum +1;
      nEvent.pfHT    += jpt;
      
      // Keep central jets
      if ( !(fabs(jeta)<3.0 ) ) continue;
      counters[PJ].fill("... |eta|<3.0");
      
      // Flag good jets failing ID
      if (!isJetID) { 
        nEvent.badJet = 1;
      } else {
        counters[PJ].fill("... pass Jet ID");
      }
      nEvent.pfGoodHT += jpt;
      sumOfPFJets += aJet;
      
      lepton tmpLepton;
      tmpLepton.p = aJet;
      tmpLepton.charge = 0;
      tmpLepton.index = i;
      tmpLepton.type = -1;
      pfGoodJets.push_back(tmpLepton);
      
      if ( jpt>30 ) {
        counters[PJ].fill("... pass tight jet selection");
        nEvent.pfTightHT += jpt;
        nEvent.pfJetGoodPt[nEvent.pfJetGoodNum]  = jpt;
        nEvent.pfJetGoodEta[nEvent.pfJetGoodNum] = jeta;
        nEvent.pfJetGoodPhi[nEvent.pfJetGoodNum] = jphi;
        nEvent.pfJetGoodID[nEvent.pfJetGoodNum]  = isJetID;
	nEvent.bTagProbCSVBP[nEvent.pfJetGoodNum] = fTR->JnewPFCombinedSecondaryVertexBPFJetTags[i];
	nEvent.bTagProbCSVMVA[nEvent.pfJetGoodNum] = fTR->JnewPFCombinedSecondaryVertexMVABPFJetTags[i];
        
        if(isJetID>0) {
	  nEvent.pfJetGoodNumID++;
	  if(nEvent.bTagProbCSVBP[nEvent.pfJetGoodNum] > 0.679) nEvent.pfJetGoodNumIDBtag++;
	}
        nEvent.pfJetGoodNum++;
	if(nEvent.bTagProbCSVBP[nEvent.pfJetGoodNum] > 0.679) nEvent.pfJetGoodNumBtag++;
      }
      if ( jpt*(jesC+unc)/jesC>30 )  nEvent.pfJetGoodNump1sigma++;
      if ( jpt*(jesC-unc)/jesC>30 )  nEvent.pfJetGoodNumn1sigma++;

      if ( jpt>40. )  nEvent.pfJetGoodNum40++;
      if ( jpt>50. )  nEvent.pfJetGoodNum50++;
      if ( jpt>60. )  nEvent.pfJetGoodNum60++;
    }
    
  TLorentzVector s1 = sortedGoodLeptons[PosLepton1].p;
  TLorentzVector s2 = sortedGoodLeptons[PosLepton2].p;

  nEvent.met[RAW]=fTR->RawMET;
  nEvent.met[DUM]=0.; // Not there anymore: fTR->MuJESCorrMET;
  nEvent.met[TCMET]=fTR->PFType1MET;
  nEvent.met[MUJESCORRMET]=fTR->MuJESCorrMET;
  nEvent.met[PFMET]=fTR->PFMET;
  nEvent.met[SUMET]=fTR->SumEt;

  TLorentzVector pfJetVector(0,0,0,0); // for constructing SumJPt from pf jets, as Pablo
  TLorentzVector pfNoCutsJetVector(0,0,0,0); // for constructing SumJPt from pfmet (unclustered), as Kostas
  TLorentzVector type1NoCutsJetVector(0,0,0,0); // for constructing SumJPt from tcmet (unclustered), new
  nEvent.metPhi[RAW]=0.;//kicked! caloMETvector.Phi();
  nEvent.metPhi[DUM]=0.;
  nEvent.metPhi[TCMET]=type1METvector.Phi();
  nEvent.metPhi[MUJESCORRMET]=0.;
  nEvent.metPhi[PFMET]=pfMETvector.Phi();
  nEvent.metPhi[SUMET]=0.;

  // remove the leptons from PFMET and tcMET blublu
  pfNoCutsJetVector = -pfMETvector - s1 - s2;
  type1NoCutsJetVector = -type1METvector - s1 - s2;

  // #--- different versions of JZB
  nEvent.dphi_sumJetVSZ[CALOJZB]=0.; // kicked caloVector.DeltaPhi(s1+s2); // DPhi between Z and SumJpt
  nEvent.sumJetPt[CALOJZB]=0.; // kicked caloVector.Pt();
  nEvent.jzb[CALOJZB] = 0.; //kicked caloVector.Pt() - (s1+s2).Pt(); // calib issue of rawcalomet wrt lepton energy scale, under develop
    
  nEvent.dphi_sumJetVSZ[PFJZB] = pfNoCutsJetVector.DeltaPhi(s1+s2); 
  nEvent.sumJetPt[PFJZB] = pfNoCutsJetVector.Pt(); 
  nEvent.jzb[PFJZB] = pfNoCutsJetVector.Pt() - (s1+s2).Pt(); // to be used with pfMET
  nEvent.sjzb[PFJZB] = GausRandom(nEvent.jzb[1]+1.3,7); // to be used with pfMET

  nEvent.dphi_sumJetVSZ[RECOILJZB] = 0.; // kicked recoil.DeltaPhi(s1+s2);  // recoil is not yet a recoil but the sumJPt, since the leptons will be added only later (ugly)
  nEvent.sumJetPt[RECOILJZB] = 0.;//kicked recoil.Pt(); 
  nEvent.jzb[RECOILJZB] = 0.;//kicked recoil.Pt() - (s1+s2).Pt(); // to be used recoil met (recoilpt[0])    
  nEvent.jzb[PFRECOILJZB] = sumOfPFJets.Pt() - (s1+s2).Pt(); // to be used recoil met (recoilpt[0])
  nEvent.sumJetPt[PFRECOILJZB] = sumOfPFJets.Pt();

  nEvent.dphi_sumJetVSZ[TCJZB] = type1NoCutsJetVector.DeltaPhi(s1+s2); // tcJZB
  nEvent.sumJetPt[TCJZB] = type1NoCutsJetVector.Pt(); 
  nEvent.jzb[TCJZB] = type1NoCutsJetVector.Pt() - (s1+s2).Pt(); // to be used with tcMET

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
  if ( nEvent.pfJetGoodNum>= 2 ) {
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

  nEvent.dphiZpfMet = (s1+s2).DeltaPhi(pfMETvector);
  nEvent.dphiZs1 = (s1+s2).DeltaPhi(s1);
  nEvent.dphiZs2 = (s1+s2).DeltaPhi(s2);
  nEvent.dphiMet1 = sortedGoodLeptons[PosLepton1].p.DeltaPhi(pfMETvector);
  nEvent.dphiMet2 = sortedGoodLeptons[PosLepton2].p.DeltaPhi(pfMETvector);
  nEvent.dphitcMet1 = sortedGoodLeptons[PosLepton1].p.DeltaPhi(type1METvector);
  nEvent.dphitcMet2 = sortedGoodLeptons[PosLepton2].p.DeltaPhi(type1METvector);
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
    int i3,i4,i5;

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
  FullTree->Write();

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

const bool JZBAnalysis::IsCustomMu2012(const int index){

  // Basic muon cleaning and ID

  // Acceptance cuts
  if (!(fTR->MuPt[index] > 10) )       return false;
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
  if ( !(fTR->MuNMatches[index] > 1) )      return false;
  counters[MU].fill(" ... nMatches > 1");
  if ( !(fTR->MuNSiLayers[index] > 5) )      return false;
  counters[MU].fill(" ... nLayers > 5");


  // Vertex compatibility
  if ( !(fabs(fTR->MuD0PV[index]) < 0.02) ) return false; //still open
  counters[MU].fill(" ... D0(pv) < 0.02");
  //HPA recommendation not POG
  if ( !(fabs(fTR->MuDzPV[index]) < 0.1 ) ) return false; //still open
  counters[MU].fill(" ... DZ(pv) < 0.1");


  // Flat isolation below 20 GeV (only for synch.: we cut at 20...)
  double Iso = MuPFIso(index);
  if ( !(Iso < 0.1) ) return false;
  counters[MU].fill(" ... Iso < 0.1");


  return true;
}

const float JZBAnalysis::GetLeptonWeight(int id1, float phi1, float eta1, int id2, float phi2, float eta2, float &EffErr) {
    //this function will become more sophisticated in the future (eta & phi based efficiency)
    if(id1==id2&&id1==0) {
      //ee
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
    }
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
  
  if(!(fabs(fTR->ElEta[index]) < 2.4) ) return false;
  counters[EL].fill(" ... |eta| < 2.4");

  if(!(fabs(fTR->ElPt[index]) > 10.0 ) ) return false;
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
  if(!(abs(fTR->ElDzPV[index])<0.1)) return false;
  counters[EL].fill(" ... DZ(PV)<0.1");

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
  double minPt = 20.;
  double mllCut = 20.;
  double maxEta = 2.4;
  double minJPt = 30;
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
  for ( int jIndex=0; jIndex<fTR->NGenJets; ++jIndex) {
    if ( fTR->GenJetPt[jIndex]<minJPt ) continue;
    if ( fabs(fTR->GenJetEta[jIndex])>maxJEta ) continue;
    ++nEvent.genNjets;
    if ( fabs(fTR->GenJetEta[jIndex])>3.0 ) continue;
    ++nEvent.genNjetsTwoSix;
  }

  // Select the two highest-pt leptons compatible with a Z
  size_t i1 = 0, i2 = 0;
  if ( sortedGLeptons.size()>1 ) {
    for ( size_t i=0; i<sortedGLeptons.size()-1; ++i ) {
      i1 = i;
      TLorentzVector lp1(sortedGLeptons[i].p);
      for ( size_t j=i+1; j<sortedGLeptons.size(); ++j ) {
        i2 = j;
        TLorentzVector lp2(sortedGLeptons[j].p);
        if ( fabs( (lp1+lp2).M() - 91.2 ) < mllCut ) break;
      }
    }
  }
  

  if(sortedGLeptons.size()>0)
    {
      nEvent.genPt1     = sortedGLeptons[i1].p.Pt();
      nEvent.genId1     = sortedGLeptons[i1].type;
      nEvent.genEta1    = sortedGLeptons[i1].p.Eta();
      nEvent.genMID     = fTR->GenLeptonMID[sortedGLeptons[i1].index];
      nEvent.genGMID    = fTR->GenLeptonGMID[sortedGLeptons[i1].index];
      if(sortedGLeptons.size()>1)
        {
          TLorentzVector genZvector = sortedGLeptons[i1].p + sortedGLeptons[i2].p;
          nEvent.genRecoil  = (-GenMETvector - genZvector).Pt();
          nEvent.genPt2     = sortedGLeptons[i2].p.Pt();
          nEvent.genId2     = sortedGLeptons[i2].type;
          nEvent.genEta2    = sortedGLeptons[i2].p.Eta();
          nEvent.genZPt     = genZvector.Pt();
          nEvent.genMll     = genZvector.M();
          nEvent.genJZB     = nEvent.genRecoil - genZvector.Pt();
	  nEvent.dphigenZgenMet = (sortedGLeptons[i1].p + sortedGLeptons[i2].p).DeltaPhi(GenMETvector);
        }
    }
}
