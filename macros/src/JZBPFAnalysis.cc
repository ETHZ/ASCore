#include "helper/Utilities.hh"
#include "JZBPFAnalysis.hh"
#include "TF1.h"
#include <time.h>
#include <TRandom.h>
#include "TLorentzVector.h"
#include "TF1.h"
//#include "/shome/theofil/setTDRStyle.C"

using namespace std;


#define jMax 30  // do not touch this
#define metMax 30
#define rMax 30

const int particleflowtypes=3+1;//  this is pf1,pf2,pf3 -- all of them get saved.  (the +1 is so that we can access pf1 with pfX[1] instead of [0] 

string sjzbPFversion="$Revision: 1.11 $";
string sjzbPFinfo="";

/*

$Log: JZBPFAnalysis.cc,v $
Revision 1.11  2011/06/09 16:28:31  buchmann
Removed jet-lepton cleaning for PF

Revision 1.10  2011/06/09 15:55:07  buchmann
Checked second occurrence. All fine

Revision 1.9  2011/06/09 15:49:41  buchmann
Fixed a bug causing Full Cleaning to lead to a segfault

Revision 1.8  2011/06/09 12:45:35  buchmann
Now also handling events that contain PF leptons but no RECO leptons

Revision 1.7  2011/06/09 12:17:36  buchmann
Adapted the getRecoElIndex and getRecoMuIndex methods to return the best candidate (it previously returned -1 if no candidate was found within 0.05, leading to a possible problem)

Revision 1.6  2011/06/09 10:44:58  pablom

Only the information of reco leptons matched to pf leptons is stored.

Revision 1.5  2011/06/09 08:54:26  buchmann
Translated more abstract variables to PF as well

Revision 1.4  2011/06/09 07:52:40  buchmann
paranthesis fixed

Revision 1.3  2011/06/09 07:46:47  buchmann
Adapted cuts and added counters for PF

Revision 1.2  2011/06/08 16:18:13  buchmann
Merged the two JZBs such that only one file is produced

Revision 1.1  2011/06/08 15:20:49  buchmann
First commit of the JZB PF analysis scripts. ATM they are practically identical to the main ones but that will change really soon.


*/


Double_t GaussRandom(Double_t mu, Double_t sigma) { 
  return gRandom->Gaus(mu,sigma);   //real deal
  //return mu;//debugging : no smearing
}
class nanoPFEvent
{
public:
  nanoPFEvent();
  void reset();

  float mll; // di-lepton system
  float recomll;
  float pfmll[particleflowtypes];
  float pt;
  float recopt;
  float pfpt[particleflowtypes];
  float phi;
  float recophi;
  float eta;
  float recoeta;
  float pfphi[particleflowtypes];
  float pfeta[particleflowtypes];
  bool is_data;

  float pt1; // leading leptons
  float pt2;
  float pfpt1[particleflowtypes]; // leading leptons
  float pfpt2[particleflowtypes];

  float recopt1; // leading leptons
  float recopt2;

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
  float recoeta1;
  float pfeta1[particleflowtypes];
  float eta2;
  float recoeta2;
  float pfeta2[particleflowtypes];
  float phi1;
  float recophi1;
  float pfphi1[particleflowtypes];
  float phi2;
  float recophi2;
  float pfphi2[particleflowtypes];
  float dphi;
  float dphiZpfMet;
  float dphiZs1;
  float dphiZs2;
  float dphiMet1;
  float dphiMet2;
  float dphitcMet1;
  float dphitcMet2;
  float dphipfRecoilMet1;
  float dphipfRecoilMet2;
  int id1;
  int id2;
  int ch1;
  int ch2;
  int chid1; // old id (kostas convention)
  int chid2;

  int jetNum;
  int goodJetNum;
  float jetpt[jMax]; // jets in barrel + endcaps
  float jeteta[jMax];
  float jetphi[jMax];
  float jetscale[jMax];
  int jetID[jMax];


  int leptonNum; // store all leptons (reduntant for the 2 leptons that form the Z)
  float leptonPt[jMax]; 
  float leptonEta[jMax];
  float leptonPhi[jMax];
  int leptonId[jMax];
  int leptonCharge[jMax];

  int pfLeptonNum; // store all leptons (reduntant for the 2 leptons that form the Z)
  float pfLeptonPt[jMax]; 
  float pfLeptonEta[jMax];
  float pfLeptonPhi[jMax];
  int pfLeptonId[jMax];
  int pfLeptonCharge[jMax];

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
  float pfJetDphiMet[jMax];
  float pfHT;
  float pfGoodHT;
  float pfTightHT;

  int pfJetGoodNum;
  int pfJetGoodNumID;
  int pfJetGoodNumEta2p4;
  int pfJetGoodNumEta2p0;
  int pfJetGoodNumEta1p4;
  int pfJetGoodNumEta1p2;
  float pfJetGoodPt[jMax];
  float pfJetGoodEta[jMax];
  float pfJetGoodPhi[jMax];
  bool   pfJetGoodID[jMax];

  int pfJetGoodNum35;
  int pfJetGoodNum33;
  int pfJetGoodNum315;
  int pfJetGoodNum285;
  int pfJetGoodNum27;
  int pfJetGoodNum25;
  int pfJetGoodNum20;


  float recoilpt[rMax];
  float dphiRecoilLep[rMax];
  float vjetpt[rMax];
  float vjeteta[rMax];
  float vjetphi[rMax];
  float recoilenergy[rMax];
  float recoilphi[rMax];
  float recoileta[rMax];

  float met[metMax];
  float metPhi[metMax];
  float dphiMetLep[metMax];
  float dphiMetJet[metMax];
  float dphiMetSumJetPt[metMax];
  float metPerp[metMax];
  float metPar[metMax];
  int eventNum;
  int runNum;
  int lumi;
  int goodVtx;
  int numVtx;
  float totEvents; // tot events processed by the ntuple producer (job submission efficiency), no need to keep this as int, better as float
  int badJet;

  float jzb[rMax];
  float pfjzb[rMax];
  float sjzb[rMax]; // smeared JZB
  float dphi_sumJetVSZ[rMax];
  float sumJetPt[rMax];

  float weight;
  float PUweight;
  bool passed_triggers;
  int trigger_bit;

};

nanoPFEvent::nanoPFEvent(){};
void nanoPFEvent::reset()
{

  mll=0; // di-lepton system
  recomll=0; // di-lepton system
  pt=0;
  recopt=0;
  phi=0;
  eta=0;
  recophi=0;

  is_data=false;

  pt1=0;
  recopt1=0;
  pt2=0;
  recopt2=0;
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
  recoeta1=0; // leading leptons
  eta2=0;
  recoeta2=0;
  phi1=0;
  recophi1=0;
  phi2=0;
  recophi2=0;
  dphiZpfMet=0;
  dphiZs1=0;
  dphiZs2=0;
  dphiMet1=0;
  dphiMet2=0;
  dphitcMet1=0;
  dphitcMet2=0;
  dphipfRecoilMet1=0;
  dphipfRecoilMet2=0;
  dphi=0;
  id1=0;
  id2=0;
  ch1=0;
  ch2=0;
  chid1=0;
  chid2=0;

  for(int ipf=0;ipf<particleflowtypes;ipf++) {
    pfmll[ipf]=0;
    pfpt[ipf]=0;
    pfpt1[ipf]=0;
    pfpt2[ipf]=0;
    pfphi[ipf]=0;
    pfphi1[ipf]=0;
    pfphi2[ipf]=0;
    pfeta[ipf]=0;
    pfeta1[ipf]=0;
    pfeta2[ipf]=0;
  }

  for(int jCounter=0;jCounter<jMax;jCounter++){
    jetpt[jCounter]=0; // jets in barrel + endcaps
    jeteta[jCounter]=0;
    jetphi[jCounter]=0;
    jetscale[jCounter]=0;
    jetID[jCounter]=0;
  }
  jetNum=0;
  goodJetNum=0;

  for(int jCounter=0;jCounter<jMax;jCounter++){
    leptonPt[jCounter]=0; 
    leptonEta[jCounter]=0;
    leptonPhi[jCounter]=0;
    leptonId[jCounter]=0;
    leptonCharge[jCounter]=0;
  }
  leptonNum=0;

  for(int jCounter=0;jCounter<jMax;jCounter++){
    pfLeptonPt[jCounter]=0; 
    pfLeptonEta[jCounter]=0;
    pfLeptonPhi[jCounter]=0;
    pfLeptonId[jCounter]=0;
    pfLeptonCharge[jCounter]=0;
  }
  pfLeptonNum=0;

  for(int jCounter=0;jCounter<jMax;jCounter++){
    leptonPairMass[jCounter]=0;
    leptonPairDphi[jCounter]=0;
    leptonPairId[jCounter]=0;
  } 
  leptonPairNum=0;
 
  for(int rCounter=0;rCounter<rMax;rCounter++){
    recoilpt[rCounter]=0;
    dphiRecoilLep[rCounter]=0;
    vjetpt[rCounter]=0;
    vjeteta[rCounter]=0;
    vjetphi[rCounter]=0;
    recoilenergy[rCounter]=0;
    recoilphi[rCounter]=0;
    recoileta[rCounter]=0;
  }

    
  for(int metCounter=0;metCounter<metMax;metCounter++){
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
  }
  pfJetGoodNum=0;
  pfJetGoodNumID=0;
  pfJetGoodNumEta2p4=0;
  pfJetGoodNumEta2p0=0;
  pfJetGoodNumEta1p4=0;
  pfJetGoodNumEta1p2=0;

  pfJetGoodNum20=0;
  pfJetGoodNum25=0;
  pfJetGoodNum27=0;
  pfJetGoodNum285=0;
  pfJetGoodNum315=0;
  pfJetGoodNum33=0;
  pfJetGoodNum35=0;

  eventNum=0;
  runNum=0;
  lumi=0;
  goodVtx=0;
  numVtx=0;
  badJet=0;
  totEvents=0;

  for(int rCounter=0;rCounter<rMax;rCounter++){
    jzb[rCounter]=0;
    sjzb[rCounter]=0;
    dphi_sumJetVSZ[rCounter]=0;
    sumJetPt[rCounter]=0;
  }

  weight = 1.0;
  PUweight = 1.0;
}


TTree *mypfTree;
TTree *InfoTree;

nanoPFEvent npfEvent;


JZBPFAnalysis::JZBPFAnalysis(TreeReader *tr, std::string dataType, bool fullCleaning) : 
  UserAnalysisBase(tr), fDataType_(dataType), fFullCleaning_(fullCleaning) {
  //	Util::SetStyle();	
  //	setTDRStyle();	
}

JZBPFAnalysis::~JZBPFAnalysis(){
}

void JZBPFAnalysis::Begin(TFile *f){
cout << endl << endl;
cout << "This is JZBPFAnalysis --- For PF: We need to decide on what a Custom PF El/Mu is (cuts); We need to define whether we want to reject events with less than 2 pf leptons or just drop back to the current way of doing things in that case (inconsistent!). We also need to look at how we want to mix the current and the new situation as with reco we may reject events that would pass pf and the other way around ... " << endl;
cout << endl << endl;
  // Define the output file of histograms
  //fHistFile = new TFile(outputFileName_.c_str(), "RECREATE");
  fHistFile=f;	
  rand_ = new TRandom();
  TH1::AddDirectory(kFALSE);

  // Define the histograms
  fHElectronPtEta = new TH2F("fHElectronPtEta","fHElectronPtEta",300,0,300,300,-3.0,3.0);
  fHElectronIDPtEta  = new TH2F("fHElectronIDPtEta","fHElectronIDPtEta",300,0,300,300,-3.0,3.0);
  fHElectronIDIsoPtEta  = new TH2F("fHElectronIDIsoPtEta","fHElectronIDIsoPtEta",300,0,300,300,-3.0,3.0);
  for(int i =0;i<20;i++){string title = "fHMee"+any2string(i);fHMee[i] = new TH1F(title.c_str(),title.c_str(),200,0,200);}
  fHMeeDPhi  = new TH2F("fHMeeDPhi","fHMeeDPhi",300,0,300,350,0,3.5);
  fHMeePt    = new TH2F("fHMeePt","fHMeePt",300,0,300,300,0,100);;
  fHMDPhiPt  = new TH2F("fHMDPhiPt","fHMDPhiPt",350,0,3.5,300,0,100);
  fHMZPtJ1Pt = new TH2F("fHMZPtJ1Pt","fHMZPtJ1Pt",300,0,300,300,0,300);
  fHMZPtuJ1Pt = new TH2F("fHMZPtuJ1Pt","fHMZPtuJ1Pt",300,0,300,300,0,300);

  InfoTree = new TTree("PFinfo","PFinfo/S");
  TString *user = new TString();
  TString *timestamp = new TString();
  TString *jzbversion = new TString();
  TString *cmsdir = new TString();
  TString *jzbinfo = new TString();
  InfoTree->Branch("user",&user,16000,0);
  InfoTree->Branch("timestamp",&timestamp,16000,0);
  InfoTree->Branch("version",&jzbversion,16000,0);
  InfoTree->Branch("cmsdir",&cmsdir,16000,0);
  InfoTree->Branch("jzbinfo",&jzbinfo,16000,0);
  char usertext[255];
  FILE *usernamefile;
  usernamefile = popen("whoami", "r");
  fgets(usertext, sizeof(usertext), usernamefile);
  pclose(usernamefile);
  *jzbversion=sjzbPFversion;
  char scmsdir[1000];
  getcwd(scmsdir,1000);
  *cmsdir=scmsdir;
  *jzbinfo=sjzbPFinfo;
  *user=usertext;
  time_t rawtime;
  struct tm * timeinfo;
  time (&rawtime );
  *timestamp=ctime(&rawtime);
  
  InfoTree->Fill();
  InfoTree->Write();

  mypfTree = new TTree("PFevents","PFevents");

  mypfTree->Branch("is_data",&npfEvent.is_data,"is_data/B");
  mypfTree->Branch("mll",&npfEvent.mll,"mll/F");
  mypfTree->Branch("recomll",&npfEvent.recomll,"recomll/F");
  mypfTree->Branch("pt",&npfEvent.pt,"pt/F");
  mypfTree->Branch("recopt",&npfEvent.recopt,"recopt/F");
  mypfTree->Branch("phi",&npfEvent.phi,"phi/F");
  mypfTree->Branch("recophi",&npfEvent.recophi,"recophi/F");
  mypfTree->Branch("pt1",&npfEvent.pt1,"pt1/F");
  mypfTree->Branch("recopt1",&npfEvent.recopt1,"recopt1/F");
  mypfTree->Branch("pt2",&npfEvent.pt2,"pt2/F");
  mypfTree->Branch("recopt2",&npfEvent.recopt2,"recopt2/F");
  mypfTree->Branch("genPt1",&npfEvent.genPt1,"genPt1/F");
  mypfTree->Branch("genPt2",&npfEvent.genPt2,"genPt2/F");
  mypfTree->Branch("genEta1",&npfEvent.genEta1,"genEta1/F");
  mypfTree->Branch("genEta2",&npfEvent.genEta2,"genEta2/F");
  mypfTree->Branch("genId1",&npfEvent.genId1,"genId1/I");
  mypfTree->Branch("genId2",&npfEvent.genId2,"genId2/I");
  mypfTree->Branch("genMID",&npfEvent.genMID,"genMID/I");
  mypfTree->Branch("genMID1",&npfEvent.genMID1,"genMID1/I");
  mypfTree->Branch("genMID2",&npfEvent.genMID2,"genMID2/I");
  mypfTree->Branch("genMID3",&npfEvent.genMID3,"genMID3/I");
  mypfTree->Branch("genMID4",&npfEvent.genMID4,"genMID4/I");
  mypfTree->Branch("genMID5",&npfEvent.genMID5,"genMID5/I");
  mypfTree->Branch("genGMID",&npfEvent.genGMID,"genGMID/I");
  mypfTree->Branch("genGMID1",&npfEvent.genGMID1,"genGMID1/I");
  mypfTree->Branch("genGMID2",&npfEvent.genGMID2,"genGMID2/I");
  mypfTree->Branch("genMET",&npfEvent.genMET,"genMET/F");
  mypfTree->Branch("genZPt",&npfEvent.genZPt,"genZPt/F");
  mypfTree->Branch("genMll",&npfEvent.genMll,"genMll/F");
  mypfTree->Branch("genRecoil",&npfEvent.genRecoil,"genRecoil/F");
  mypfTree->Branch("genJZB",&npfEvent.genJZB,"genJZB/F");
  mypfTree->Branch("genNjets",&npfEvent.genNjets,"genNjets/I");
  mypfTree->Branch("genNleptons",&npfEvent.genNleptons,"genNleptons/I");
  mypfTree->Branch("genPt1Sel",&npfEvent.genPt1Sel,"genPt1Sel/F");
  mypfTree->Branch("genPt2Sel",&npfEvent.genPt2Sel,"genPt2Sel/F");
  mypfTree->Branch("genEta1Sel",&npfEvent.genEta1Sel,"genEta1Sel/F");
  mypfTree->Branch("genEta2Sel",&npfEvent.genEta2Sel,"genEta2Sel/F");
  mypfTree->Branch("genId1Sel",&npfEvent.genId1Sel,"genId1Sel/I");
  mypfTree->Branch("genId2Sel",&npfEvent.genId2Sel,"genId2Sel/I");
  mypfTree->Branch("genZPtSel",&npfEvent.genZPtSel,"genZPtSel/F");
  mypfTree->Branch("genMllSel",&npfEvent.genMllSel,"genMllSel/F");
  mypfTree->Branch("genRecoilSel",&npfEvent.genRecoilSel,"genRecoilSel/F");
  mypfTree->Branch("genJZBSel",&npfEvent.genJZBSel,"genJZBSel/F");
  mypfTree->Branch("eta1",&npfEvent.eta1,"eta1/F");
  mypfTree->Branch("recoeta1",&npfEvent.recoeta1,"recoeta1/F");
  mypfTree->Branch("eta2",&npfEvent.eta2,"eta2/F");
  mypfTree->Branch("recoeta2",&npfEvent.recoeta2,"recoeta2/F");
  mypfTree->Branch("phi1",&npfEvent.phi1,"phi1/F");
  mypfTree->Branch("recophi1",&npfEvent.recophi1,"recophi1/F");
  mypfTree->Branch("phi2",&npfEvent.phi2,"phi2/F");
  mypfTree->Branch("recophi2",&npfEvent.recophi2,"recophi2/F");
  mypfTree->Branch("dphiZpfMet",&npfEvent.dphiZpfMet,"dphiZpfMet/F");
  mypfTree->Branch("dphiZs1",&npfEvent.dphiZs1,"dphiZs1/F");
  mypfTree->Branch("dphiZs2",&npfEvent.dphiZs2,"dphiZs2/F");
  mypfTree->Branch("dphiMet1",&npfEvent.dphiMet1,"dphiMet1/F");
  mypfTree->Branch("dphiMet2",&npfEvent.dphiMet2,"dphiMet2/F");
  mypfTree->Branch("dphitcMet1",&npfEvent.dphitcMet1,"dphitcMet1/F");
  mypfTree->Branch("dphitcMet2",&npfEvent.dphitcMet2,"dphitcMet2/F");
  mypfTree->Branch("dphipfRecoilMet1",&npfEvent.dphipfRecoilMet1,"dphipfRecoilMet1/F");
  mypfTree->Branch("dphipfRecoilMet2",&npfEvent.dphipfRecoilMet2,"dphipfRecoilMet2/F");
  mypfTree->Branch("dphi",&npfEvent.dphi,"dphi/F");

  mypfTree->Branch("id1",&npfEvent.id1,"id1/I");
  mypfTree->Branch("id2",&npfEvent.id2,"id2/I");
  mypfTree->Branch("ch1",&npfEvent.ch1,"ch1/I");
  mypfTree->Branch("ch2",&npfEvent.ch2,"ch2/I");
  mypfTree->Branch("chid1",&npfEvent.chid1,"chid1/I");
  mypfTree->Branch("chid2",&npfEvent.chid2,"chid2/I");

  mypfTree->Branch("jetNum",&npfEvent.jetNum,"jetNum/I");
  mypfTree->Branch("goodJetNum",&npfEvent.goodJetNum,"goodJetNum/I");
  mypfTree->Branch("jetID",npfEvent.jetID,"jetID[jetNum]/I");
  mypfTree->Branch("jetpt",npfEvent.jetpt,"jetpt[jetNum]/F");
  mypfTree->Branch("jeteta",npfEvent.jeteta,"jeteta[jetNum]/F");
  mypfTree->Branch("jetphi",npfEvent.jetphi,"jetphi[jetNum]/F");
  mypfTree->Branch("jetscale",npfEvent.jetscale,"jetscale[jetNum]/F");

  mypfTree->Branch("leptonNum",&npfEvent.leptonNum,"leptonNum/I");
  mypfTree->Branch("leptonPt",npfEvent.leptonPt,"leptonPt[leptonNum]/F");
  mypfTree->Branch("leptonEta",npfEvent.leptonEta,"leptonEta[leptonNum]/F");
  mypfTree->Branch("leptonPhi",npfEvent.leptonPhi,"leptonPhi[leptonNum]/F");
  mypfTree->Branch("leptonId",npfEvent.leptonId,"leptonId[leptonNum]/I");
  mypfTree->Branch("leptonCharge",npfEvent.leptonCharge,"leptonCharge[leptonNum]/I");

  mypfTree->Branch("pfLeptonNum",&npfEvent.pfLeptonNum,"pfLeptonNum/I");
  mypfTree->Branch("pfLeptonPt",npfEvent.pfLeptonPt,"pfLeptonPt[pfLeptonNum]/F");
  mypfTree->Branch("pfLeptonEta",npfEvent.pfLeptonEta,"pfLeptonEta[pfLeptonNum]/F");
  mypfTree->Branch("pfLeptonPhi",npfEvent.pfLeptonPhi,"pfLeptonPhi[pfLeptonNum]/F");
  mypfTree->Branch("pfLeptonId",npfEvent.pfLeptonId,"pfLeptonId[pfLeptonNum]/I");
  mypfTree->Branch("pfLeptonCharge",npfEvent.pfLeptonCharge,"pfLeptonCharge[pfLeptonNum]/I");

  mypfTree->Branch("leptonPairNum",&npfEvent.leptonPairNum,"leptonPairNum/I");
  mypfTree->Branch("leptonPairMass",npfEvent.leptonPairMass,"leptonPairMass[leptonPairNum]/F");
  mypfTree->Branch("leptonPairDphi",npfEvent.leptonPairDphi,"leptonPairDphi[leptonPairNum]/F");
  mypfTree->Branch("leptonPairId",npfEvent.leptonPairId,"leptonPairId[leptonPairNum]/I");

  mypfTree->Branch("recoilpt",npfEvent.recoilpt,"recoilpt[30]/F");
  mypfTree->Branch("dphiRecoilLep",npfEvent.dphiRecoilLep,"dphiRecoilLep[30]/F");
  mypfTree->Branch("recoilphi",npfEvent.recoilphi,"recoilphi[30]/F");
  mypfTree->Branch("recoileta",npfEvent.recoileta,"recoileta[30]/F");
  mypfTree->Branch("recoilenergy",npfEvent.recoilenergy,"recoilenergy[30]/F");

  mypfTree->Branch("vjetpt",npfEvent.vjetpt,"vjetpt[30]/F");
  mypfTree->Branch("vjeteta",npfEvent.vjeteta,"vjeteta[30]/F");
  mypfTree->Branch("vjetphi",npfEvent.vjetphi,"vjetphi[30]/F");

  mypfTree->Branch("met",npfEvent.met,"met[30]/F");
  mypfTree->Branch("metPhi",npfEvent.metPhi,"metPhi[30]/F");
  mypfTree->Branch("dphiMetLep",npfEvent.dphiMetLep,"dphiMetLep[30]/F");
  mypfTree->Branch("dphiMetJet",npfEvent.dphiMetJet,"dphiMetJet[30]/F");
  mypfTree->Branch("dphiMetSumJetPt",npfEvent.dphiMetSumJetPt,"dphiMetSumJetPt[30]/F");
  mypfTree->Branch("metPerp",npfEvent.metPerp,"metPerp[30]/F");
  mypfTree->Branch("metPar",npfEvent.metPar,"metPar[30]/F");


  mypfTree->Branch("eventNum",&npfEvent.eventNum,"eventNum/I");
  mypfTree->Branch("runNum",&npfEvent.runNum,"runNum/I");
  mypfTree->Branch("lumi",&npfEvent.lumi,"lumi/I");
  mypfTree->Branch("goodVtx",&npfEvent.goodVtx,"goodVtx/I");
  mypfTree->Branch("numVtx",&npfEvent.numVtx,"numVtx/I");
  mypfTree->Branch("badJet",&npfEvent.badJet,"badJet/I");
  mypfTree->Branch("totEvents",&npfEvent.totEvents,"totEvents/F");

  mypfTree->Branch("pfJetNum",&npfEvent.pfJetNum,"pfJetNum/I");
  mypfTree->Branch("pfJetPt",npfEvent.pfJetPt,"pfJetPt[pfJetNum]/F");
  mypfTree->Branch("pfJetEta",npfEvent.pfJetEta,"pfJetEta[pfJetNum]/F");
  mypfTree->Branch("pfJetPhi",npfEvent.pfJetPhi,"pfJetPhi[pfJetNum]/F");
  mypfTree->Branch("pfJetID",npfEvent.pfJetID,"pfJetID[pfJetNum]/B");
  mypfTree->Branch("pfJetScale",npfEvent.pfJetScale,"pfJetScale[pfJetNum]/F");
  mypfTree->Branch("pfJetDphiMet",npfEvent.pfJetDphiMet,"pfJetDphiMet[pfJetNum]/F");
  mypfTree->Branch("pfHT",&npfEvent.pfHT,"pfHT/F");
  mypfTree->Branch("pfGoodHT",&npfEvent.pfGoodHT,"pfGoodHT/F");
  mypfTree->Branch("pfTightHT",&npfEvent.pfTightHT,"pfTightHT/F");

  mypfTree->Branch("pfJetGoodNum",&npfEvent.pfJetGoodNum,"pfJetGoodNum/I");
  mypfTree->Branch("pfJetGoodNumID",&npfEvent.pfJetGoodNumID,"pfJetGoodNumID/I");
  mypfTree->Branch("pfJetGoodNumEta2p4",&npfEvent.pfJetGoodNumEta2p4,"pfJetGoodNumEta2p4/I");
  mypfTree->Branch("pfJetGoodNumEta2p0",&npfEvent.pfJetGoodNumEta2p0,"pfJetGoodNumEta2p0/I");
  mypfTree->Branch("pfJetGoodNumEta1p4",&npfEvent.pfJetGoodNumEta1p4,"pfJetGoodNumEta1p4/I");
  mypfTree->Branch("pfJetGoodNumEta1p2",&npfEvent.pfJetGoodNumEta1p2,"pfJetGoodNumEta1p2/I");

  mypfTree->Branch("pfJetGoodPt", npfEvent.pfJetGoodPt,"pfJetGoodPt[pfJetGoodNum]/F");
  mypfTree->Branch("pfJetGoodEta",npfEvent.pfJetGoodEta,"pfJetGoodEta[pfJetGoodNum]/F");
  mypfTree->Branch("pfJetGoodPhi",npfEvent.pfJetGoodPhi,"pfJetGoodPhi[pfJetGoodNum]/F");
  mypfTree->Branch("pfJetGoodID", npfEvent.pfJetGoodID,"pfJetGoodID[pfJetGoodNum]/B");

  mypfTree->Branch("pfJetGoodNum20",&npfEvent.pfJetGoodNum20,"pfJetGoodNum20/I");
  mypfTree->Branch("pfJetGoodNum25",&npfEvent.pfJetGoodNum25,"pfJetGoodNum25/I");
  mypfTree->Branch("pfJetGoodNum27",&npfEvent.pfJetGoodNum27,"pfJetGoodNum27/I");
  mypfTree->Branch("pfJetGoodNum285",&npfEvent.pfJetGoodNum285,"pfJetGoodNum285/I");
  mypfTree->Branch("pfJetGoodNum315",&npfEvent.pfJetGoodNum315,"pfJetGoodNum315/I");
  mypfTree->Branch("pfJetGoodNum33",&npfEvent.pfJetGoodNum33,"pfJetGoodNum33/I");
  mypfTree->Branch("pfJetGoodNum35",&npfEvent.pfJetGoodNum35,"pfJetGoodNum35/I");

  mypfTree->Branch("pfpt",&npfEvent.pfpt,"pfpt[30]/F");
  mypfTree->Branch("pfphi",&npfEvent.pfphi,"pfphi[30]/F");
  mypfTree->Branch("pfeta",&npfEvent.pfeta,"pfeta[30]/F");
  mypfTree->Branch("eta",&npfEvent.eta,"eta/F");
  mypfTree->Branch("recoeta",&npfEvent.recoeta,"recoeta/F");
  mypfTree->Branch("recoeta",&npfEvent.recoeta,"recoeta/F");
  mypfTree->Branch("pfpt1",&npfEvent.pfpt1,"pfpt1[30]/F");
  mypfTree->Branch("pfpt2",&npfEvent.pfpt2,"pfpt2[30]/F");
  mypfTree->Branch("pfphi1",&npfEvent.pfphi1,"pfphi1[30]/F");
  mypfTree->Branch("pfphi2",&npfEvent.pfphi2,"pfphi2[30]/F");
  mypfTree->Branch("pfeta1",&npfEvent.pfeta1,"pfeta1[30]/F");
  mypfTree->Branch("pfeta2",&npfEvent.pfeta2,"pfeta2[30]/F");
  mypfTree->Branch("pfmll",&npfEvent.pfmll,"pfmll[30]/F");
  mypfTree->Branch("pfjzb",&npfEvent.pfjzb,"pfjzb[30]/F");

  mypfTree->Branch("jzb",npfEvent.jzb,"jzb[30]/F");
  mypfTree->Branch("sjzb",npfEvent.sjzb,"sjzb[30]/F");
  mypfTree->Branch("dphi_sumJetVSZ",npfEvent.dphi_sumJetVSZ,"dphi_sumJetVSZ[30]/F");
  mypfTree->Branch("sumJetPt",npfEvent.sumJetPt,"sumJetPt[30]/F");
  mypfTree->Branch("weight", &npfEvent.weight,"weight/F");
  mypfTree->Branch("PUweight",&npfEvent.PUweight,"PUweight/F");
  mypfTree->Branch("passed_triggers", &npfEvent.passed_triggers,"passed_triggers/B");
  mypfTree->Branch("trigger_bit", &npfEvent.trigger_bit,"trigger_bit/I");


  // Define counters (so we have them in the right order)
  counters[EV].setName("Events");
  counters[TR].setName("Triggers");
  counters[MU].setName("Muons");
  counters[PFMU].setName("PFMuons");
  counters[EL].setName("Electrons");
  counters[PFEL].setName("PFElectrons");
  counters[JE].setName("Jets");
  counters[PJ].setName("PFJets");

  counters[EV].fill("All events",0.);
  if ( fDataType_ != "mc" ) {
    counters[EV].fill("... pass electron triggers",0.);
    counters[EV].fill("... pass muon triggers",0.);
    counters[EV].fill("... pass EM triggers",0.);
    counters[EV].fill("... pass all trigger requirements",0.);
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
bool momentumComparator(PFlepton i, PFlepton j) { return (i.p.Pt()>j.p.Pt()); }


//------------------------------------------------------------------------------
vector<PFlepton> JZBPFAnalysis::sortLeptonsByPt(vector<PFlepton>& leptons) {
  
  vector<PFlepton> theLep = leptons;
  sort (theLep.begin(), theLep.end(), momentumComparator);
  return theLep;  
  
}


//------------------------------------------------------------------------------
const bool JZBPFAnalysis::passElTriggers() {
  if ( GetHLTResult("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v1") )        return true;
  if ( GetHLTResult("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v2") )        return true;
  if ( GetHLTResult("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v3") )        return true;
  if ( GetHLTResult("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v4") )        return true;
  if ( GetHLTResult("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v5") )        return true;
  return false;

}

//------------------------------------------------------------------------------
const bool JZBPFAnalysis::passMuTriggers() {
  if ( GetHLTResult("HLT_DoubleMu6_v1") )        return true;
  if ( GetHLTResult("HLT_DoubleMu6_v2") )        return true;
  if ( GetHLTResult("HLT_DoubleMu6_v3") )        return true;
  if ( GetHLTResult("HLT_DoubleMu7_v1") )        return true;
  if ( GetHLTResult("HLT_DoubleMu7_v2") )        return true;
  if ( GetHLTResult("HLT_DoubleMu7_v3") )        return true;
  return false;
} 

//______________________________________________________________________________
const bool JZBPFAnalysis::passEMuTriggers() {
  if ( GetHLTResult("HLT_Mu17_Ele8_CaloIdL_v1") )        return true;
  if ( GetHLTResult("HLT_Mu17_Ele8_CaloIdL_v2") )        return true;
  if ( GetHLTResult("HLT_Mu17_Ele8_CaloIdL_v3") )        return true;
  if ( GetHLTResult("HLT_Mu17_Ele8_CaloIdL_v4") )        return true;
  if ( GetHLTResult("HLT_Mu17_Ele8_CaloIdL_v5") )        return true;
  if ( GetHLTResult("HLT_Mu8_Ele17_CaloIdL_v1") )        return true;
  if ( GetHLTResult("HLT_Mu8_Ele17_CaloIdL_v2") )        return true;
  if ( GetHLTResult("HLT_Mu8_Ele17_CaloIdL_v3") )        return true;
  if ( GetHLTResult("HLT_Mu8_Ele17_CaloIdL_v4") )        return true;
  if ( GetHLTResult("HLT_Mu8_Ele17_CaloIdL_v5") )        return true;
  return false;
}


//______________________________________________________________________________
void JZBPFAnalysis::Analyze() {

  // #--- analysis global parameters
  double DRmax=0.4; // veto jets in a cone of DRmax close to the lepton

  counters[EV].fill("All events");
  npfEvent.reset();
  // Fill generic information
  npfEvent.eventNum  = fTR->Event;
  npfEvent.runNum    = fTR->Run;
  npfEvent.lumi      = fTR->LumiSection;
  npfEvent.totEvents = fTR->GetEntries();
  if(fDataType_ == "mc") // only do this for MC; for data npfEvent.reset() has already set both weights to 1 
    {
      npfEvent.PUweight  = GetPUWeight(fTR->PUnumInteractions);
      npfEvent.weight    = GetPUWeight(fTR->PUnumInteractions);
    }
  // Trigger information
  npfEvent.passed_triggers=0;
  if ( fDataType_ != "mc" ) 
    {
      npfEvent.is_data=true;
      if ( passElTriggers() ) 
        {
          counters[EV].fill("... pass electron triggers");
          npfEvent.passed_triggers=1;
          npfEvent.trigger_bit |= 1;
        } 
      if ( passMuTriggers() )
        {
          counters[EV].fill("... pass muon triggers");
          npfEvent.passed_triggers=1;
          npfEvent.trigger_bit |= (1<<1);
        } 
      if ( passEMuTriggers() )
        {
          counters[EV].fill("... pass EM triggers");
          npfEvent.passed_triggers=1;
          npfEvent.trigger_bit |= (1<<2);
        }
    }

  // Check if we find an OSSF pair in the acceptance (and if it is coming from a Z)
  bool isMC = (fDataType_ == "mc");
  if ( isMC ) { GeneratorInfo(); }
    
  // #--- Vertex info
  npfEvent.numVtx = fTR->NVrtx;
  float rho = sqrt(fTR->PrimVtxx*fTR->PrimVtxx + fTR->PrimVtxy*fTR->PrimVtxy);
  if(fTR->PrimVtxGood) npfEvent.goodVtx |=2; // save bits of vertex quality
  if (   fTR->PrimVtxGood==0 && fTR->PrimVtxIsFake==0 
      && fTR->PrimVtxNdof>4  && fabs(fTR->PrimVtxz)<24 && rho<2)
    npfEvent.goodVtx |=4;
  
  // Good event requirement: essentially vertex requirements
  if ( !IsGoodEvent() ) {
    if( isMC ) mypfTree->Fill();
    return;
  }
  counters[EV].fill("... pass good event requirements");

  vector<PFlepton> leptons;
  vector<vector<PFlepton> > pfLeptons;

  for(int ipf=0;ipf<particleflowtypes;ipf++) {
	vector<PFlepton> pfLepton;
	pfLeptons.push_back(pfLepton);
  }

  TLorentzVector genZvector; // To store the true Z vector

//--------------------------   STEP 1 : RECO LEPTONS

  // #--- muon loop
  /*  for(int muIndex=0;muIndex<fTR->NMus;muIndex++)
    {
      counters[MU].fill("All mus");
      if(IsCustomMu(muIndex))
	{
	  counters[MU].fill("... pass mu selection");
	  float px= fTR->MuPx[muIndex];
	  float py= fTR->MuPy[muIndex];
	  float pz= fTR->MuPz[muIndex];
	  float energy =  fTR->MuE[muIndex];
	  TLorentzVector tmpVector(px,py,pz,energy);
	  int tmpCharge = fTR->MuCharge[muIndex];
	  
	  PFlepton tmpLepton;
	  tmpLepton.p = tmpVector;
	  tmpLepton.charge = tmpCharge;
	  tmpLepton.index = muIndex;
	  tmpLepton.type = 1;
	  tmpLepton.genPt = 0.;
	  leptons.push_back(tmpLepton);
	}
    }
  
  
  // #--- electron loop
  for(int elIndex=0;elIndex<fTR->NEles;elIndex++)
    {
      counters[EL].fill("All eles");
      if(IsCustomEl(elIndex))	
	{
          counters[EL].fill("... pass e selection");
	  float px= fTR->ElPx[elIndex];
	  float py= fTR->ElPy[elIndex];
	  float pz= fTR->ElPz[elIndex];
	  float energy =  fTR->ElE[elIndex];
	  TLorentzVector tmpVector(px,py,pz,energy);
	  int tmpCharge=fTR->ElCharge[elIndex];
	  
	  PFlepton tmpLepton;
	  tmpLepton.p = tmpVector;
	  tmpLepton.charge = tmpCharge;
	  tmpLepton.index = elIndex;
	  tmpLepton.type = 0;
	  tmpLepton.genPt = 0.;
	  leptons.push_back(tmpLepton);
	}
    }
   */

//--------------------------   STEP 2 : PF LEPTONS

  // #-- PF muon loop -- type 1
  for(int muIndex=0;muIndex<fTR->PfMuNObjs;muIndex++)
    {
      if(IsCustomPfMu(muIndex)) {
	  TLorentzVector tmpVector(fTR->PfMuPx[muIndex],fTR->PfMuPy[muIndex],fTR->PfMuPz[muIndex],fTR->PfMuE[muIndex]);
	  int recoIndex = getRecoMuIndex(tmpVector); //Get the index of a matched reco muon
	  TLorentzVector recotmpVector(0,0,0,0);
	  int recotmpCharge = 0;
	  if(recoIndex>-1) {
		recotmpVector = (fTR->MuPx[recoIndex],fTR->MuPy[recoIndex],fTR->MuPz[recoIndex],fTR->MuE[recoIndex]);
		recotmpCharge = fTR->MuCharge[recoIndex];
	  }
	  int tmpCharge = fTR->PfMu3Charge[muIndex];
	
          PFlepton tmpLepton;
	  tmpLepton.p = tmpVector;
          tmpLepton.recop = recotmpVector;
	  tmpLepton.charge = tmpCharge;
	  tmpLepton.recocharge = recotmpCharge;
	  tmpLepton.index = muIndex;
	  tmpLepton.recoindex = recoIndex;
	  tmpLepton.type = 1;
	  tmpLepton.genPt = 0.;
	  pfLeptons[1].push_back(tmpLepton);
      }
    }


  // #-- PF muon loop -- type 2
  for(int muIndex=0;muIndex<fTR->PfMu2NObjs;muIndex++)
    {
      counters[PFMU].fill("All PF mus");
      if ( npfEvent.pfLeptonNum>=jMax ) break;
      npfEvent.pfLeptonPt[npfEvent.pfLeptonNum]     = fTR->PfMu3Pt[muIndex];
      npfEvent.pfLeptonEta[npfEvent.pfLeptonNum]    = fTR->PfMu3Eta[muIndex];
      npfEvent.pfLeptonPhi[npfEvent.pfLeptonNum]    = fTR->PfMu3Phi[muIndex];
      npfEvent.pfLeptonId[npfEvent.pfLeptonNum]     = 13*fTR->PfMu2Charge[muIndex];
      npfEvent.pfLeptonCharge[npfEvent.pfLeptonNum] = fTR->PfMu2Charge[muIndex];
      npfEvent.pfLeptonNum++;
      if ( npfEvent.pfLeptonNum>=jMax ) break;
      if(IsCustomPfMu(muIndex)) {
	  TLorentzVector tmpVector(fTR->PfMu2Px[muIndex],fTR->PfMu2Py[muIndex],fTR->PfMu2Pz[muIndex],fTR->PfMu2E[muIndex]);
	  int recoIndex = getRecoMuIndex(tmpVector); //Get the index of a matched reco muon
	  TLorentzVector recotmpVector(0,0,0,0);
	  int recotmpCharge = 0;
	  if(recoIndex>-1) {
		recotmpVector = (fTR->MuPx[recoIndex],fTR->MuPy[recoIndex],fTR->MuPz[recoIndex],fTR->MuE[recoIndex]);
		recotmpCharge = fTR->MuCharge[recoIndex];
	  }
	  int tmpCharge = fTR->PfMu3Charge[muIndex];

	  PFlepton tmpLepton;
	  tmpLepton.p = tmpVector;
	  tmpLepton.recop = recotmpVector;
	  tmpLepton.charge = tmpCharge;
	  tmpLepton.recocharge = recotmpCharge;
	  tmpLepton.index = muIndex;
	  tmpLepton.recoindex = recoIndex;
	  tmpLepton.type = 1;
	  tmpLepton.genPt = 0.;
	  pfLeptons[2].push_back(tmpLepton);
	  pfLeptons[0].push_back(tmpLepton); // THIS IS THE REAL DEAL (the one we will probably want to use)
      }
    }

  // #-- PF muon loop -- type 3
  for(int muIndex=0;muIndex<fTR->PfMu3NObjs;muIndex++)
    {
      if(IsCustomPfMu(muIndex)) {
          counters[MU].fill("... pass pf mu selection");
	  TLorentzVector tmpVector(fTR->PfMu3Px[muIndex],fTR->PfMu3Py[muIndex],fTR->PfMu3Pz[muIndex],fTR->PfMu3E[muIndex]);
	  int recoIndex = getRecoMuIndex(tmpVector); //Get the index of a matched reco muon
	  TLorentzVector recotmpVector(0,0,0,0);
	  int recotmpCharge = 0;
	  if(recoIndex>-1) {
		recotmpVector = (fTR->MuPx[recoIndex],fTR->MuPy[recoIndex],fTR->MuPz[recoIndex],fTR->MuE[recoIndex]);
		recotmpCharge = fTR->MuCharge[recoIndex];
	  }
	  int tmpCharge = fTR->PfMu3Charge[muIndex];

	  PFlepton tmpLepton;
	  tmpLepton.p = tmpVector;
	  tmpLepton.recop = recotmpVector;
	  tmpLepton.charge = tmpCharge;
	  tmpLepton.recocharge = recotmpCharge;
	  tmpLepton.index = muIndex;
	  tmpLepton.recoindex = recoIndex;
	  tmpLepton.type = 1;
	  tmpLepton.genPt = 0.;
	  pfLeptons[3].push_back(tmpLepton);
      }
    }


  // #-- PF electron loop -- type 1
  for(int elIndex=0;elIndex<fTR->PfElNObjs;elIndex++)
    {
      counters[PFEL].fill("All PF eles");
      if(IsCustomPfEl(elIndex)) {
      if ( npfEvent.pfLeptonNum>=jMax ) break;
      npfEvent.pfLeptonPt[npfEvent.pfLeptonNum]     = fTR->PfElPt[elIndex];
      npfEvent.pfLeptonEta[npfEvent.pfLeptonNum]    = fTR->PfElEta[elIndex];
      npfEvent.pfLeptonPhi[npfEvent.pfLeptonNum]    = fTR->PfElPhi[elIndex];
      npfEvent.pfLeptonId[npfEvent.pfLeptonNum]     = 11*fTR->PfElCharge[elIndex];
      npfEvent.pfLeptonCharge[npfEvent.pfLeptonNum] = fTR->PfElCharge[elIndex];
      npfEvent.pfLeptonNum++;
	  TLorentzVector tmpVector(fTR->PfElPx[elIndex],fTR->PfElPy[elIndex],fTR->PfElPz[elIndex],fTR->PfElE[elIndex]);
	  int recoIndex = getRecoElIndex(tmpVector); //Get the index of a matched reco electron
	  TLorentzVector recotmpVector(0,0,0,0);
	  int recotmpCharge = 0;
	  if(recoIndex>-1) {
		recotmpVector = (fTR->ElPx[recoIndex],fTR->ElPy[recoIndex],fTR->ElPz[recoIndex],fTR->ElE[recoIndex]);
		recotmpCharge = fTR->ElCharge[recoIndex];
	  }
	  int tmpCharge = fTR->PfMu3Charge[elIndex];

	  PFlepton tmpLepton;
	  tmpLepton.p = tmpVector;
	  tmpLepton.recop = recotmpVector;
	  tmpLepton.charge = tmpCharge;
	  tmpLepton.recocharge = recotmpCharge;
	  tmpLepton.index = elIndex;
	  tmpLepton.recoindex = recoIndex;
	  tmpLepton.type = 0;
	  tmpLepton.genPt = 0.;
	  pfLeptons[1].push_back(tmpLepton);
	  pfLeptons[0].push_back(tmpLepton); // THIS IS THE REAL DEAL (the one we will probably want to use)
      }
    }
  
  // #-- PF electron loop -- type 2
  for(int elIndex=0;elIndex<fTR->PfEl2NObjs;elIndex++)
    {
      if(IsCustomPfEl(elIndex)) {
	  TLorentzVector tmpVector(fTR->PfEl2Px[elIndex],fTR->PfEl2Py[elIndex],fTR->PfEl2Pz[elIndex],fTR->PfEl2E[elIndex]);
	  int recoIndex = getRecoElIndex(tmpVector); //Get the index of a matched reco electron
	  TLorentzVector recotmpVector(0,0,0,0);
	  int recotmpCharge = 0;
	  if(recoIndex>-1) {
		recotmpVector = (fTR->ElPx[recoIndex],fTR->ElPy[recoIndex],fTR->ElPz[recoIndex],fTR->ElE[recoIndex]);
		recotmpCharge = fTR->ElCharge[recoIndex];
	  }
	  int tmpCharge = fTR->PfEl3Charge[elIndex];
	 
	  PFlepton tmpLepton;
	  tmpLepton.p = tmpVector;
	  tmpLepton.recop = recotmpVector;
	  tmpLepton.charge = tmpCharge;
	  tmpLepton.recocharge = recotmpCharge;
	  tmpLepton.index = elIndex;
	  tmpLepton.recoindex = recoIndex;
	  tmpLepton.type = 0;
	  tmpLepton.genPt = 0.;
	  pfLeptons[2].push_back(tmpLepton);
      }
    }
  
  // #-- PF electron loop -- type 3
  for(int elIndex=0;elIndex<fTR->PfEl3NObjs;elIndex++)
    {
      if(IsCustomPfEl(elIndex)) {
          counters[EL].fill("... pass pf e selection");
	  TLorentzVector tmpVector(fTR->PfEl3Px[elIndex],fTR->PfEl3Py[elIndex],fTR->PfEl3Pz[elIndex],fTR->PfEl3E[elIndex]);
	  int recoIndex = getRecoElIndex(tmpVector); //Get the index of a matched reco electron
	
	  TLorentzVector recotmpVector(0,0,0,0);
	  int recotmpCharge = 0;
	  if(recoIndex>-1) {
		recotmpVector = (fTR->ElPx[recoIndex],fTR->ElPy[recoIndex],fTR->ElPz[recoIndex],fTR->ElE[recoIndex]);
		recotmpCharge = fTR->ElCharge[recoIndex];
	  }
	  int tmpCharge = fTR->PfEl3Charge[elIndex];

	  PFlepton tmpLepton;
	  tmpLepton.p = tmpVector;
	  tmpLepton.recop = recotmpVector;
	  tmpLepton.charge = tmpCharge;
	  tmpLepton.recocharge = recotmpCharge;
	  tmpLepton.index = elIndex;
	  tmpLepton.recoindex = recoIndex;
	  tmpLepton.type = 0;
	  tmpLepton.genPt = 0.;
	  pfLeptons[3].push_back(tmpLepton);
      }
    }

  // Sort the leptons by Pt and select the two opposite-signed ones with highest Pt

  //vector<PFlepton> sortedGoodLeptons = sortLeptonsByPt(leptons);
 
  vector<vector<PFlepton> > sortedGoodPFLeptons;
  bool dopf[particleflowtypes];
  
  //bool doreco;
  
  for(int ipf=0;ipf<particleflowtypes;ipf++) {
	vector<PFlepton> leptonsel = pfLeptons[ipf];
	vector<PFlepton> sortedGoodPfLeptonOfSpecificPFtype;
	if(leptonsel.size()>0) sortedGoodPfLeptonOfSpecificPFtype =sortLeptonsByPt(leptonsel);
        sortedGoodPFLeptons.push_back(sortedGoodPfLeptonOfSpecificPFtype);
        //if(sortedGoodPFLeptons[ipf].size() > 1) {dopf[ipf]=true;} else {dopf[ipf]=false;}
  }
  
  
  //if(sortedGoodLeptons.size() > 1) {doreco=true;} else {doreco=false;}

  if(sortedGoodPFLeptons[0].size() > 1) { // note that the "0" entry corresponds to the one we actually want to use.
    counters[EV].fill("... has at least 2 leptons");
    //int PosLepton1 = 0;
    //int PosLepton2 = 1;
    int PfPosLepton1[particleflowtypes];
    int PfPosLepton2[particleflowtypes];
    for(int ipf=0;ipf<particleflowtypes;ipf++) {
	PfPosLepton1[ipf]=0;
	PfPosLepton2[ipf]=1;
    }

    // Check for OS combination
    //for(; PosLepton2 < sortedGoodLeptons.size(); PosLepton2++) {
    //  if(sortedGoodLeptons[0].charge*sortedGoodLeptons[PosLepton2].charge<0) break;
    //}
    for(int ipf=0;ipf<particleflowtypes;ipf++) {
      for(; PfPosLepton2[ipf] < sortedGoodPFLeptons[ipf].size(); PfPosLepton2[ipf]++) {
	if(sortedGoodPFLeptons[ipf][0].charge*sortedGoodPFLeptons[ipf][PfPosLepton2[ipf]].charge<0) break;
      }
    }
	
    if(PfPosLepton2[0] == sortedGoodPFLeptons[0].size()) {//not enough leptons...
      if(isMC) mypfTree->Fill(); //only fill it if we're dealing with MC
      return;
    }
    counters[EV].fill("... has at least 2 OS leptons");
    // is this a smart way of doing it?
    
    // Preselection
    /*if(sortedGoodLeptons.size()>1 && PosLepton2!=sortedGoodLeptons.size() && sortedGoodLeptons[PosLepton1].p.Pt() > 20 && sortedGoodLeptons[PosLepton2].p.Pt() > 20) {

      npfEvent.eta1 = sortedGoodLeptons[PosLepton1].p.Eta();
      npfEvent.pt1 = sortedGoodLeptons[PosLepton1].p.Pt();
      npfEvent.phi1 = sortedGoodLeptons[PosLepton1].p.Phi();
      npfEvent.ch1 = sortedGoodLeptons[PosLepton1].charge;
      npfEvent.id1 = sortedGoodLeptons[PosLepton1].type; //??????
      npfEvent.chid1 = (sortedGoodLeptons[PosLepton1].type+1)*sortedGoodLeptons[PosLepton1].charge;
      
      npfEvent.eta2 = sortedGoodLeptons[PosLepton2].p.Eta();
      npfEvent.pt2 = sortedGoodLeptons[PosLepton2].p.Pt();
      npfEvent.phi2 = sortedGoodLeptons[PosLepton2].p.Phi();
      npfEvent.ch2 = sortedGoodLeptons[PosLepton2].charge;
      npfEvent.id2 = sortedGoodLeptons[PosLepton2].type; //??????
      npfEvent.chid2 = (sortedGoodLeptons[PosLepton2].type+1)*sortedGoodLeptons[PosLepton2].charge;

      npfEvent.mll=(sortedGoodLeptons[PosLepton2].p+sortedGoodLeptons[PosLepton1].p).M();
      npfEvent.phi=(sortedGoodLeptons[PosLepton2].p+sortedGoodLeptons[PosLepton1].p).Phi();
      npfEvent.pt=(sortedGoodLeptons[PosLepton2].p+sortedGoodLeptons[PosLepton1].p).Pt();
      npfEvent.dphi=sortedGoodLeptons[PosLepton2].p.DeltaPhi(sortedGoodLeptons[PosLepton1].p);
    } else {
      
      //If there are less than two leptons the data is not saved.
      mypfTree->Fill();
      return;
      
    }
    */

    if(!(sortedGoodPFLeptons[0][PfPosLepton1[0]].p.Pt() > 20 && sortedGoodPFLeptons[0][PfPosLepton2[0]].p.Pt() > 20)) {
      if(isMC) mypfTree->Fill();
      return;
    }
    counters[EV].fill("... pass dilepton pt selection");
    
    //Fil in PF information
    npfEvent.pt1=sortedGoodPFLeptons[0][PfPosLepton1[0]].p.Pt();
    npfEvent.pt2=sortedGoodPFLeptons[0][PfPosLepton2[0]].p.Pt();
    npfEvent.phi1=sortedGoodPFLeptons[0][PfPosLepton1[0]].p.Phi();
    npfEvent.phi2=sortedGoodPFLeptons[0][PfPosLepton2[0]].p.Phi();
    npfEvent.eta1=sortedGoodPFLeptons[0][PfPosLepton1[0]].p.Eta();
    npfEvent.eta2=sortedGoodPFLeptons[0][PfPosLepton2[0]].p.Eta();
    npfEvent.mll=(sortedGoodPFLeptons[0][PfPosLepton1[0]].p+sortedGoodPFLeptons[0][PfPosLepton2[0]].p).M();
    npfEvent.phi=(sortedGoodPFLeptons[0][PfPosLepton1[0]].p+sortedGoodPFLeptons[0][PfPosLepton2[0]].p).Phi();
    npfEvent.eta=(sortedGoodPFLeptons[0][PfPosLepton1[0]].p+sortedGoodPFLeptons[0][PfPosLepton2[0]].p).Eta();
    npfEvent.pt=(sortedGoodPFLeptons[0][PfPosLepton1[0]].p+sortedGoodPFLeptons[0][PfPosLepton2[0]].p).Pt();

    //Fill the reco information
    npfEvent.recopt1=sortedGoodPFLeptons[0][PfPosLepton1[0]].recop.Pt();
    npfEvent.recopt2=sortedGoodPFLeptons[0][PfPosLepton2[0]].recop.Pt();
    npfEvent.recophi1=sortedGoodPFLeptons[0][PfPosLepton1[0]].recop.Phi();
    npfEvent.recophi2=sortedGoodPFLeptons[0][PfPosLepton2[0]].recop.Phi();
    npfEvent.recoeta1=sortedGoodPFLeptons[0][PfPosLepton1[0]].recop.Eta();
    npfEvent.recoeta2=sortedGoodPFLeptons[0][PfPosLepton2[0]].recop.Eta();
    npfEvent.recomll=(sortedGoodPFLeptons[0][PfPosLepton1[0]].recop+sortedGoodPFLeptons[0][PfPosLepton2[0]].recop).M();
    npfEvent.recophi=(sortedGoodPFLeptons[0][PfPosLepton1[0]].recop+sortedGoodPFLeptons[0][PfPosLepton2[0]].recop).Phi();
    npfEvent.recoeta=(sortedGoodPFLeptons[0][PfPosLepton1[0]].recop+sortedGoodPFLeptons[0][PfPosLepton2[0]].recop).Eta();
    npfEvent.recopt=(sortedGoodPFLeptons[0][PfPosLepton1[0]].recop+sortedGoodPFLeptons[0][PfPosLepton2[0]].recop).Pt();

    for(int ipf=0;ipf<particleflowtypes;ipf++) {
      if(dopf[ipf]&&sortedGoodPFLeptons[ipf][PfPosLepton1[ipf]].p.Pt()>20&&sortedGoodPFLeptons[ipf][PfPosLepton2[ipf]].p.Pt()>20) {
	npfEvent.pfpt1[ipf]=sortedGoodPFLeptons[ipf][PfPosLepton1[ipf]].p.Pt();
	npfEvent.pfpt2[ipf]=sortedGoodPFLeptons[ipf][PfPosLepton2[ipf]].p.Pt();
	npfEvent.pfphi1[ipf]=sortedGoodPFLeptons[ipf][PfPosLepton1[ipf]].p.Phi();
	npfEvent.pfphi2[ipf]=sortedGoodPFLeptons[ipf][PfPosLepton2[ipf]].p.Phi();
	npfEvent.pfeta1[ipf]=sortedGoodPFLeptons[ipf][PfPosLepton1[ipf]].p.Eta();
	npfEvent.pfeta2[ipf]=sortedGoodPFLeptons[ipf][PfPosLepton2[ipf]].p.Eta();
	npfEvent.pfmll[ipf]=(sortedGoodPFLeptons[ipf][PfPosLepton1[ipf]].p+sortedGoodPFLeptons[ipf][PfPosLepton2[ipf]].p).M();
	npfEvent.pfphi[ipf]=(sortedGoodPFLeptons[ipf][PfPosLepton1[ipf]].p+sortedGoodPFLeptons[ipf][PfPosLepton2[ipf]].p).Phi();
	npfEvent.pfphi[ipf]=(sortedGoodPFLeptons[ipf][PfPosLepton1[ipf]].p+sortedGoodPFLeptons[ipf][PfPosLepton2[ipf]].p).Eta();
	npfEvent.pfpt[ipf]=(sortedGoodPFLeptons[ipf][PfPosLepton1[ipf]].p+sortedGoodPFLeptons[ipf][PfPosLepton2[ipf]].p).Pt();
      }
    }
    
    
    // #--- construct different recoil models, initial the recoil vector will hold only the sum over the hard jets, only in the end we will add-up the lepton system
    TLorentzVector recoil(0,0,0,0); // different constructions of recoil model (under dev, need cleaning)    
    npfEvent.jetNum=0;        // total jet counting
    npfEvent.goodJetNum=0;    // Jets passing tighter pt cut
    for(int i =0 ; i<fTR->NJets;i++) // CALO jet loop
      {
        counters[JE].fill("All Calo jets");
	if(i==jMax) { cout<<"max Num was reached"<<endl; return; }
	
	float jpt  = fTR->JPt[i];
	float jeta = fTR->JEta[i];
	float jpx  = fTR->JPx[i];
	float jpy  = fTR->JPy[i];
	float jpz  = fTR->JPz[i];
	float jenergy = fTR->JE[i];
	float jesC    = fTR->JEcorr[i];
	bool isJetID  = IsCustomJet(i);
	
        // Consider only Jets passing JetID
	if (!isJetID) continue;
        //FIXME: throw away event if good jet does not pass JetID?
        counters[JE].fill("... pass jet ID");
	
	TLorentzVector aJet(jpx,jpy,jpz,jenergy);

/*
//Superfluous now ... was only there to avoid duplication when using both PF jets and reco leptons.
        // lepton-jet cleaning
        if ( fFullCleaning_ ) { 
          // Remove jet close to any lepton
          bool isClean(true);
          for ( size_t ilep = 0; ilep<sortedGoodPFLeptons[0].size(); ++ilep )
            if ( aJet.DeltaR(sortedGoodPFLeptons[0][ilep].p)<DRmax) isClean=false;
          if ( !isClean ) continue;
          counters[JE].fill("... pass full lepton cleaning");
        } else {
          // Remove jet close to leptons from Z candidate
          if(aJet.DeltaR(sortedGoodPFLeptons[0][PfPosLepton1[0]].p)<DRmax)continue; 
          counters[JE].fill("... pass lepton 1 veto");
          if(aJet.DeltaR(sortedGoodPFLeptons[0][PfPosLepton1[0]].p)<DRmax)continue;
          counters[JE].fill("... pass lepton 2 veto");
        }
*/	
        // Acceptance cuts before we use this jet
        if ( !(fabs(jeta)<2.6 && jpt>20.) ) continue;
        counters[JE].fill("... |eta|<2.6 && pt>20.");
        recoil+=aJet;
	
	npfEvent.jetpt[npfEvent.jetNum]  = aJet.Pt();
	npfEvent.jeteta[npfEvent.jetNum] = aJet.Eta();
	npfEvent.jetphi[npfEvent.jetNum] = aJet.Phi();
	npfEvent.jetscale[npfEvent.jetNum]  = jesC;
	if(isJetID) npfEvent.jetID[npfEvent.jetNum] = 1;
	
	npfEvent.jetNum = npfEvent.jetNum + 1 ;
	
        if ( jpt>30 ) {
          counters[JE].fill("... pt>30");
          npfEvent.goodJetNum++;
        }
	
      }
    

    // --- construct met vectors here
    float caloMETpx = fTR->RawMETpx;
    float caloMETpy = fTR->RawMETpy;
    
    float pfMETpx = fTR->PFMETpx;
    float pfMETpy = fTR->PFMETpy;

    float tcMETpx = fTR->TCMETpx;
    float tcMETpy = fTR->TCMETpy;
    
    TLorentzVector caloMETvector(caloMETpx,caloMETpy,0,0);
    TLorentzVector pfMETvector(pfMETpx,pfMETpy,0,0);
    TLorentzVector tcMETvector(tcMETpx,tcMETpy,0,0);
    TLorentzVector sumOfPFJets(0,0,0,0);

    vector<PFlepton> pfGoodJets;
    for(int i =0 ; i<fTR->NJets;i++) // PF jet loop//killPF
      {
        counters[PJ].fill("All PF jets");
	if(i==jMax){cout<<"max Num was reached"<<endl; return;}
	
	float jpt = fTR->JPt[i];//killPF
	float jeta = fTR->JEta[i];//killPF
	float jphi = fTR->JPhi[i];//killPF
	float jpx = fTR->JPx[i];//killPF
	float jpy = fTR->JPy[i];//killPF
	float jpz = fTR->JPz[i];//killPF
	float jenergy = fTR->JE[i];//killPF
	float jesC = fTR->JEcorr[i];//killPF
	bool  isJetID = IsGoodBasicPFJet(i,false);
	
	TLorentzVector aJet(jpx,jpy,jpz,jenergy);

/*
	//Superfluous now ... was only there to avoid duplication when using both PF jets and reco leptons.
        // lepton-jet cleaning
        if ( fFullCleaning_ ) { 
          // Remove jet close to any lepton
          bool isClean(true);
          for ( size_t ilep = 0; ilep<sortedGoodPFLeptons[0].size(); ++ilep )
            if ( aJet.DeltaR(sortedGoodPFLeptons[0][ilep].p)<DRmax) isClean=false;
          if ( !isClean ) continue;
          counters[PJ].fill("... pass full lepton cleaning");
        } else {
          // Remove jet close to leptons from Z candidate
          if(aJet.DeltaR(sortedGoodPFLeptons[0][PfPosLepton1[0]].p)<DRmax)continue; 
          counters[PJ].fill("... pass lepton 1 veto");
          if(aJet.DeltaR(sortedGoodPFLeptons[0][PfPosLepton2[0]].p)<DRmax)continue;
          counters[PJ].fill("... pass lepton 2 veto");
        }
*/
	
        // Keep jets over min. pt threshold
        if ( !(jpt>20) ) continue;  
        counters[PJ].fill("... pt>20.");

	npfEvent.pfJetPt[npfEvent.pfJetNum]    = jpt;
	npfEvent.pfJetEta[npfEvent.pfJetNum]   = jeta;
	npfEvent.pfJetPhi[npfEvent.pfJetNum]   = jphi;
        npfEvent.pfJetScale[npfEvent.pfJetNum] = jesC;
        npfEvent.pfJetID[npfEvent.pfJetNum]    = isJetID;
	npfEvent.pfJetDphiMet[npfEvent.pfJetNum] = aJet.DeltaPhi(pfMETvector);
	npfEvent.pfJetNum = npfEvent.pfJetNum +1;
	npfEvent.pfHT    += jpt;

        // Keep central jets
        if ( !(fabs(jeta)<2.6 ) ) continue;
        counters[PJ].fill("... |eta|<2.6");

        // Flag good jets failing ID
	if (!isJetID) { 
          npfEvent.badJet = 1;
        } else {
          counters[PJ].fill("... pass Jet ID");
        }
        npfEvent.pfGoodHT += jpt;
        sumOfPFJets += aJet;

        PFlepton tmpLepton;
        tmpLepton.p = aJet;
        tmpLepton.charge = 0;
        tmpLepton.index = i;
        tmpLepton.type = -1;
        pfGoodJets.push_back(tmpLepton);

        if ( jpt>30 ) {
          counters[PJ].fill("... pass tight jet selection");
          npfEvent.pfTightHT += jpt;
          npfEvent.pfJetGoodPt[npfEvent.pfJetGoodNum]  = jpt;
          npfEvent.pfJetGoodEta[npfEvent.pfJetGoodNum] = jeta;
          npfEvent.pfJetGoodPhi[npfEvent.pfJetGoodNum] = jphi;
          npfEvent.pfJetGoodID[npfEvent.pfJetGoodNum]  = isJetID;
          if(isJetID>0) npfEvent.pfJetGoodNumID++;
          npfEvent.pfJetGoodNum++;
          if (abs(jeta)<2.4) npfEvent.pfJetGoodNumEta2p4++;
          if (abs(jeta)<2.0) npfEvent.pfJetGoodNumEta2p0++;
          if (abs(jeta)<1.4) npfEvent.pfJetGoodNumEta1p4++;
          if (abs(jeta)<1.2) npfEvent.pfJetGoodNumEta1p2++;
        }
        if ( jpt>20. )  npfEvent.pfJetGoodNum20++;
        if ( jpt>25. )  npfEvent.pfJetGoodNum25++;
        if ( jpt>27. )  npfEvent.pfJetGoodNum27++;
        if ( jpt>28.5 ) npfEvent.pfJetGoodNum285++;
        if ( jpt>31.5 ) npfEvent.pfJetGoodNum315++;
        if ( jpt>33. )  npfEvent.pfJetGoodNum33++;
        if ( jpt>35. )  npfEvent.pfJetGoodNum35++;
      }
    
    
    int index;
    if(recoil.Pt()!=0) // so far we had not added the lepton system in the recoil, so our recoil represents the sumJPt (ugly but it should work)
      {
	index=0;
	npfEvent.vjetpt[index]=recoil.Pt();  // vjet = vector sum of jets, vjetpt = sumJPt
	npfEvent.vjeteta[index]=recoil.Eta();
	npfEvent.vjetphi[index]=recoil.Phi();
      }
    
    /*
    deliberately commented out at the moment to make sure these two are not used anymore

    TLorentzVector s1;
    if(doreco) s1 = sortedGoodLeptons[PosLepton1].p;
    TLorentzVector s2;
    if(doreco) s2 = sortedGoodLeptons[PosLepton2].p;
*/
    TLorentzVector PFs1[particleflowtypes];
    TLorentzVector PFs2[particleflowtypes];

    for(int ipf=0;ipf<particleflowtypes;ipf++) {
	if(dopf[ipf]) PFs1[ipf] = sortedGoodPFLeptons[ipf][PfPosLepton1[ipf]].p;
	if(dopf[ipf]) PFs2[ipf] = sortedGoodPFLeptons[ipf][PfPosLepton2[ipf]].p;
    }

    npfEvent.met[0]=fTR->RawMET;
    npfEvent.met[1]=0.; // Not there anymore: fTR->MuJESCorrMET;
    npfEvent.met[2]=fTR->TCMET;
    npfEvent.met[3]=fTR->MuJESCorrMET;
    npfEvent.met[4]=fTR->PFMET;
    npfEvent.met[5]=fTR->SumEt;
//    npfEvent.met[6]=fTR->METPAT; ??

    TLorentzVector caloVector(0,0,0,0); // for constructing SumJPt from raw calomet
    TLorentzVector pfJetVector(0,0,0,0); // for constructing SumJPt from pf jets, as Pablo
    TLorentzVector pfrecoNoCutsJetVector(0,0,0,0); // for constructing SumJPt from pfmet (unclustered), as Kostas
    TLorentzVector purepfNoCutsJetVector(0,0,0,0); // for constructing SumJPt from pfmet (unclustered), as Kostas
    TLorentzVector AllpfNoCutsJetVector[particleflowtypes]; // for constructing SumJPt from pfmet (unclustered), as Kostas
    TLorentzVector tcNoCutsJetVector(0,0,0,0); // for constructing SumJPt from tcmet (unclustered), new
    npfEvent.metPhi[0]=caloMETvector.Phi();
    npfEvent.metPhi[1]=0.;
    npfEvent.metPhi[2]=tcMETvector.Phi();
    npfEvent.metPhi[3]=0.;
    npfEvent.metPhi[4]=pfMETvector.Phi();
    npfEvent.metPhi[5]=0.;
    
    
    // Remove electrons from MET
    caloVector = -caloMETvector;
//    if ( doreco && sortedGoodLeptons[PosLepton1].type == 0 ) caloVector -= s1;
//    if ( doreco && sortedGoodLeptons[PosLepton2].type == 0 ) caloVector -= s2;
    if ( sortedGoodPFLeptons[0][PfPosLepton1[0]].type == 0 ) caloVector -= PFs1[0];
    if ( sortedGoodPFLeptons[0][PfPosLepton2[0]].type == 0 ) caloVector -= PFs2[0];

    // remove the leptons from PFMET and tcMET
//    if(doreco) pfrecoNoCutsJetVector = -pfMETvector - s1 - s2;
//    if(doreco) tcNoCutsJetVector = -tcMETvector - s1 - s2;
    purepfNoCutsJetVector  = -pfMETvector - PFs1[0] - PFs2[0];
    tcNoCutsJetVector = -tcMETvector - PFs1[0] - PFs2[0];

    for(int ipf=0;ipf<particleflowtypes;ipf++) {
	if(dopf[ipf]) AllpfNoCutsJetVector[ipf] = -pfMETvector - PFs1[ipf] - PFs2[ipf];
    }
    // #--- different versions of JZB
    
    //these do not require RECO
//    npfEvent.dphi_sumJetVSZ[2] = recoil.DeltaPhi(s1+s2);  // recoil is not yet a recoil but the sumJPt, since the leptons will be added only later (ugly)
    npfEvent.dphi_sumJetVSZ[2] = recoil.DeltaPhi(PFs1[0]+PFs2[0]);  // recoil is not yet a recoil but the sumJPt, since the leptons will be added only later (ugly)
    npfEvent.sumJetPt[2] = recoil.Pt(); 
    npfEvent.sumJetPt[3] = sumOfPFJets.Pt();
    
    for(int ipf=0;ipf<particleflowtypes;ipf++) {
	if(dopf[ipf]) npfEvent.pfjzb[ipf]=AllpfNoCutsJetVector[ipf].Pt() - (PFs1[ipf]+PFs2[ipf]).Pt();
	else npfEvent.pfjzb[ipf]=-99999.9;
    }
    
    //these do require RECO
//    if(doreco) {
//      npfEvent.dphi_sumJetVSZ[0]=caloVector.DeltaPhi(s1+s2); // DPhi between Z and SumJpt
//      npfEvent.sumJetPt[0]=caloVector.Pt();
//      npfEvent.jzb[0] = caloVector.Pt() - (s1+s2).Pt(); // calib issue of rawcalomet wrt lepton energy scale, under develop
      npfEvent.dphi_sumJetVSZ[0]=caloVector.DeltaPhi(PFs1[0]+PFs2[0]); // DPhi between Z and SumJpt
      npfEvent.sumJetPt[0]=caloVector.Pt();
      npfEvent.jzb[0] = caloVector.Pt() - (PFs1[0]+PFs2[0]).Pt(); // calib issue of rawcalomet wrt lepton energy scale, under develop
      
//      npfEvent.dphi_sumJetVSZ[1] = pfrecoNoCutsJetVector.DeltaPhi(s1+s2); 
//      npfEvent.sumJetPt[1] = pfrecoNoCutsJetVector.Pt(); 
//      npfEvent.jzb[1] = pfrecoNoCutsJetVector.Pt() - (s1+s2).Pt(); // to be used with pfMET
      npfEvent.dphi_sumJetVSZ[1] = purepfNoCutsJetVector.DeltaPhi(PFs1[0]+PFs2[0]); 
      npfEvent.sumJetPt[1] = purepfNoCutsJetVector.Pt(); 
      npfEvent.jzb[1] = purepfNoCutsJetVector.Pt() - (PFs1[0]+PFs2[0]).Pt(); // to be used with pfMET

      npfEvent.sjzb[1] = GaussRandom(npfEvent.jzb[1]+1.3,7);
      
      
//      npfEvent.jzb[2] = recoil.Pt() - (s1+s2).Pt(); // to be used recoil met (recoilpt[0])    
//      npfEvent.jzb[3] = sumOfPFJets.Pt() - (s1+s2).Pt(); // to be used recoil met (recoilpt[0])
      npfEvent.jzb[2] = recoil.Pt() - (PFs1[0]+PFs2[0]).Pt(); // to be used recoil met (recoilpt[0])    
      npfEvent.jzb[3] = sumOfPFJets.Pt() - (PFs1[0]+PFs2[0]).Pt(); // to be used recoil met (recoilpt[0])
      
//      npfEvent.dphi_sumJetVSZ[4] = tcNoCutsJetVector.DeltaPhi(s1+s2); // tcJZB
//      npfEvent.sumJetPt[4] = tcNoCutsJetVector.Pt(); 
//      npfEvent.jzb[4] = tcNoCutsJetVector.Pt() - (s1+s2).Pt(); // to be used with tcMET
      npfEvent.dphi_sumJetVSZ[4] = tcNoCutsJetVector.DeltaPhi(PFs1[0]+PFs2[0]); // tcJZB
      npfEvent.sumJetPt[4] = tcNoCutsJetVector.Pt(); 
      npfEvent.jzb[4] = tcNoCutsJetVector.Pt() - (PFs1[0]+PFs2[0]).Pt(); // to be used with tcMET
      
      // --- recoil met and pf recoil met
//      npfEvent.met[6] = (sumOfPFJets + s1 + s2).Pt(); 
//      npfEvent.met[7] = (recoil + s1 + s2).Pt();
      npfEvent.met[6] = (sumOfPFJets + PFs1[0]+PFs2[0]).Pt(); 
      npfEvent.met[7] = (recoil + PFs1[0]+PFs2[0]).Pt();
      
//      recoil+=s1+s2;   // now add also the leptons to the recoil! to form the complete recoil model
      recoil+=PFs1[0]+PFs2[0];   // now add also the leptons to the recoil! to form the complete recoil model

      if(recoil.Pt()!=0)
      {
	index=0;
	npfEvent.recoilpt[index]=recoil.Pt();
	npfEvent.recoileta[index]=recoil.Eta();
	npfEvent.recoilphi[index]=recoil.Phi();
      }
//    }//end of if(doreco)
    
    // ----------------------------------------
    
    // Statistics ///////////////////////////////////////
    string type("");
    switch ( (npfEvent.id1+1)*(npfEvent.id2+1) ) {
    case 1: type = "ee"; break;
    case 2: type = "em"; break;
    case 4: type = "mm"; break;
    default: type = "unknown";
    }
    counters[EV].fill("... "+type+" pairs");     
    if ( npfEvent.pfJetGoodNum>= 2 ) {
      counters[EV].fill("... "+type+" + 2 jets");
      if ( fabs(npfEvent.mll-91)<20 ) {
        counters[EV].fill("... "+type+" + 2 jets + require Z");
        if ( npfEvent.jzb[1]>50 ) {
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
    npfEvent.leptonNum = int(sortedGoodPFLeptons[0].size());
    for ( size_t i=0; i<sortedGoodPFLeptons[0].size(); ++i ) {
      TLorentzVector lp(sortedGoodPFLeptons[0][i].p);
      npfEvent.leptonPt[i] = lp.Pt();
      npfEvent.leptonEta[i] = lp.Eta();
      npfEvent.leptonPhi[i] = lp.Phi();
      npfEvent.leptonCharge[i] = sortedGoodPFLeptons[0][i].charge;
      npfEvent.leptonId[i] = sortedGoodPFLeptons[0][i].type ;


      for(size_t j=i+1; j<sortedGoodPFLeptons[0].size();j++) // store lepton pair masses
	{
          TLorentzVector lp1(sortedGoodPFLeptons[0][i].p);
          TLorentzVector lp2(sortedGoodPFLeptons[0][j].p);
          int old_id1 = (sortedGoodPFLeptons[0][i].type+1)*sortedGoodPFLeptons[0][i].charge;
          int old_id2 = (sortedGoodPFLeptons[0][j].type+1)*sortedGoodPFLeptons[0][j].charge;
          if(npfEvent.leptonPairNum<jMax)
	    {
              npfEvent.leptonPairMass[npfEvent.leptonPairNum] = (lp1+lp2).M();
              npfEvent.leptonPairDphi[npfEvent.leptonPairNum] = lp1.DeltaPhi(lp2);
              npfEvent.leptonPairId[npfEvent.leptonPairNum] = old_id1*old_id2;
              npfEvent.leptonPairNum=npfEvent.leptonPairNum+1;
	    }
	}
    }
    /*
    if(doreco) {
      npfEvent.dphiZpfMet = (s1+s2).DeltaPhi(pfMETvector);
      npfEvent.dphiZs1 = (s1+s2).DeltaPhi(s1);
      npfEvent.dphiZs2 = (s1+s2).DeltaPhi(s2);
      npfEvent.dphiMet1 = sortedGoodLeptons[PosLepton1].p.DeltaPhi(pfMETvector);
      npfEvent.dphiMet2 = sortedGoodLeptons[PosLepton2].p.DeltaPhi(pfMETvector);
      npfEvent.dphitcMet1 = sortedGoodLeptons[PosLepton1].p.DeltaPhi(tcMETvector);
      npfEvent.dphitcMet2 = sortedGoodLeptons[PosLepton2].p.DeltaPhi(tcMETvector);
      npfEvent.dphipfRecoilMet1 = sortedGoodLeptons[PosLepton1].p.DeltaPhi(-sumOfPFJets - s1 - s2); // pf recoil met
      npfEvent.dphipfRecoilMet2 = sortedGoodLeptons[PosLepton2].p.DeltaPhi(-sumOfPFJets - s1 - s2); // pf recoil met
    }
    */

      npfEvent.dphiZpfMet = (PFs1[0]+PFs2[0]).DeltaPhi(pfMETvector);
      npfEvent.dphiZs1 = (PFs1[0]+PFs2[0]).DeltaPhi(PFs1[0]);
      npfEvent.dphiZs2 = (PFs1[0]+PFs2[0]).DeltaPhi(PFs2[0]);
      npfEvent.dphiMet1 = sortedGoodPFLeptons[0][PfPosLepton1[0]].p.DeltaPhi(pfMETvector);
      npfEvent.dphiMet2 = sortedGoodPFLeptons[0][PfPosLepton2[0]].p.DeltaPhi(pfMETvector);
      npfEvent.dphitcMet1 = sortedGoodPFLeptons[0][PfPosLepton1[0]].p.DeltaPhi(tcMETvector);
      npfEvent.dphitcMet2 = sortedGoodPFLeptons[0][PfPosLepton2[0]].p.DeltaPhi(tcMETvector);
      npfEvent.dphipfRecoilMet1 = sortedGoodPFLeptons[0][PfPosLepton1[0]].p.DeltaPhi(-sumOfPFJets - PFs1[0]-PFs2[0]); // pf recoil met
      npfEvent.dphipfRecoilMet2 = sortedGoodPFLeptons[0][PfPosLepton2[0]].p.DeltaPhi(-sumOfPFJets - PFs1[0]-PFs2[0]); // pf recoil met

    // Store minimum dphi between some mets and any kind of lepton
/*
    for ( size_t i=0; i<sortedGoodLeptons.size(); ++i ) {
      TLorentzVector lp(sortedGoodLeptons[i].p);
      if ( fabs(recoil.DeltaPhi(lp))>npfEvent.dphiRecoilLep[0] )
        npfEvent.dphiRecoilLep[0] = recoil.DeltaPhi(lp);
      if ( fabs(pfMETvector.DeltaPhi(lp))>npfEvent.dphiMetLep[4] )
        npfEvent.dphiMetLep[4] = pfMETvector.DeltaPhi(lp);
      if ( fabs((sumOfPFJets + s1 + s2).DeltaPhi(lp))> npfEvent.dphiMetLep[6] )
        npfEvent.dphiMetLep[6] = (sumOfPFJets + s1 + s2).DeltaPhi(lp);
      if ( fabs((recoil + s1 + s2).DeltaPhi(lp)) > npfEvent.dphiMetLep[7] )
        npfEvent.dphiMetLep[7] = (recoil + s1 + s2).DeltaPhi(lp);
    }
*/
    for ( size_t i=0; i<sortedGoodPFLeptons[0].size(); ++i ) {
      TLorentzVector lp(sortedGoodPFLeptons[0][i].p);
      if ( fabs(recoil.DeltaPhi(lp))>npfEvent.dphiRecoilLep[0] )
        npfEvent.dphiRecoilLep[0] = recoil.DeltaPhi(lp);
      if ( fabs(pfMETvector.DeltaPhi(lp))>npfEvent.dphiMetLep[4] )
        npfEvent.dphiMetLep[4] = pfMETvector.DeltaPhi(lp);
      if ( fabs((sumOfPFJets + PFs1[0]+PFs2[0]).DeltaPhi(lp))> npfEvent.dphiMetLep[6] )
        npfEvent.dphiMetLep[6] = (sumOfPFJets + PFs1[0]+PFs2[0]).DeltaPhi(lp);
      if ( fabs((recoil + PFs1[0]+PFs2[0]).DeltaPhi(lp)) > npfEvent.dphiMetLep[7] )
        npfEvent.dphiMetLep[7] = (recoil + PFs1[0]+PFs2[0]).DeltaPhi(lp);
    }

    // Store minimum dphi between some mets and any good jet
    for ( size_t i=0; i<pfGoodJets.size(); ++i ) {
      TLorentzVector jp(pfGoodJets[i].p);
      if ( fabs(pfMETvector.DeltaPhi(jp))>npfEvent.dphiMetJet[4] )
        npfEvent.dphiMetJet[4] = pfMETvector.DeltaPhi(jp);
    }
/*
    if(doreco) {
      // Store some additional MET information
      npfEvent.dphiMetSumJetPt[4] = pfrecoNoCutsJetVector.DeltaPhi(pfMETvector);
      npfEvent.metPar[4]  = pfMETvector.Dot(s1+s2);
      npfEvent.metPerp[4] = pfMETvector.Perp((s1+s2).Vect());
    }*/
      // Store some additional MET information
      npfEvent.dphiMetSumJetPt[4] = pfrecoNoCutsJetVector.DeltaPhi(pfMETvector);
      npfEvent.metPar[4]  = pfMETvector.Dot(PFs1[0]+PFs2[0]);
      npfEvent.metPerp[4] = pfMETvector.Perp((PFs1[0]+PFs2[0]).Vect());
    
    // Store some generator information on selected leptons
    if ( isMC ) {
      TLorentzVector GenMETvector(fTR->GenMETpx,fTR->GenMETpy,0,0);
      int i1 = sortedGoodPFLeptons[0][PfPosLepton1[0]].index;
      int i2 = sortedGoodPFLeptons[0][PfPosLepton2[0]].index;
      int i3,i4,i5;

      /// REALLY NOT SURE HERE!

      TLorentzVector genLep1; 
      if ( sortedGoodPFLeptons[0][PfPosLepton1[0]].type )
        genLep1.SetPtEtaPhiE(fTR->MuGenPt[i1],fTR->MuGenEta[i1],fTR->MuGenPhi[i1],fTR->MuGenE[i1]);
      else
        genLep1.SetPtEtaPhiE(fTR->ElGenPt[i1],fTR->ElGenEta[i1],fTR->ElGenPhi[i1],fTR->ElGenE[i1]);
      TLorentzVector genLep2;
      if ( sortedGoodPFLeptons[0][PfPosLepton2[0]].type )
        genLep2.SetPtEtaPhiE(fTR->MuGenPt[i2],fTR->MuGenEta[i2],fTR->MuGenPhi[i2],fTR->MuGenE[i2]);
      else
        genLep2.SetPtEtaPhiE(fTR->ElGenPt[i2],fTR->ElGenEta[i2],fTR->ElGenPhi[i2],fTR->ElGenE[i2]);
      
      npfEvent.genRecoilSel = (-GenMETvector - genLep1 - genLep2).Pt();
      npfEvent.genZPtSel    = (genLep1 + genLep2).Pt();
      npfEvent.genMllSel    = (genLep1 + genLep2).M();

      npfEvent.genMID1     = fTR->GenLeptonMID[i1]; // WW study
      npfEvent.genMID2     = fTR->GenLeptonMID[i2]; // WW study

      if(sortedGoodPFLeptons[0].size()>=3) {i3=sortedGoodPFLeptons[0][2].index;npfEvent.genMID3 = fTR->GenLeptonMID[i3];} else npfEvent.genMID3=-999; // WW study
      if(sortedGoodPFLeptons[0].size()>=4) {i4=sortedGoodPFLeptons[0][3].index;npfEvent.genMID4 = fTR->GenLeptonMID[i4];} else npfEvent.genMID4=-999; // WW study
      if(sortedGoodPFLeptons[0].size()>=5) {i5=sortedGoodPFLeptons[0][4].index;npfEvent.genMID5 = fTR->GenLeptonMID[i5];} else npfEvent.genMID5=-999; // WW study

      npfEvent.genGMID1    = fTR->GenLeptonGMID[i1]; // WW study
      npfEvent.genGMID2    = fTR->GenLeptonGMID[i2]; // WW study

      npfEvent.genPt1Sel    = genLep1.Pt();
      npfEvent.genPt2Sel    = genLep2.Pt();
      npfEvent.genEta1Sel   = genLep1.Eta();
      npfEvent.genEta2Sel   = genLep2.Eta();
      npfEvent.genId1Sel    = (sortedGoodPFLeptons[0][PfPosLepton1[0]].type?fTR->MuGenID[i1]:fTR->ElGenID[i1]);
      npfEvent.genId2Sel    = (sortedGoodPFLeptons[0][PfPosLepton2[0]].type?fTR->MuGenID[i2]:fTR->ElGenID[i2]);
      npfEvent.genJZBSel    = npfEvent.genRecoilSel - (genLep1 + genLep2).Pt();
	
    }
    mypfTree->Fill();
    }
}

void JZBPFAnalysis::End(TFile *f){
  f->cd();	

  mypfTree->Write();

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
//  fHistFile->Close();
}

template<class T>
std::string JZBPFAnalysis::any2string(T i)
{
  std::ostringstream buffer;
  buffer << i;
  return buffer.str();
}

const bool JZBPFAnalysis::IsCustomPfMu(const int index){
//VERY TEMPORARY !!!!
  // Basic muon cleaning and ID
  // Acceptance cuts
  if ( !(fTR->PfMu3Pt[index] > 10) )       return false;
  counters[PFMU].fill(" ... PF pt > 10");
  if ( !(fabs(fTR->PfMu3Eta[index])<2.4) ) return false;
  counters[PFMU].fill(" ... PF |eta| < 2.4");

/*  // Quality cuts
  if ( !fTR->fTR->PfMu3IsGMPT[index])        return false;
  counters[MU].fill(" ... is global muon prompt tight");
  if ( !fTR->MuIsGlobalMuon[index] )  return false;
  counters[MU].fill(" ... is global muon");
  if ( !fTR->MuIsTrackerMuon[index] ) return false;
  counters[MU].fill(" ... is tracker muon");

  // Hits
  if ( !(fTR->MuNTkHits[index] >= 11) )     return false;
  counters[MU].fill(" ... nTkHits >= 11");
  if ( !(fTR->MuNPxHits[index] > 0) )       return false;
  counters[MU].fill(" ... nPxHits > 0");
  if ( !(fTR->MuNMatches[index] > 1) )      return false;
  counters[MU].fill(" ... nMatches > 1");

  // Vertex compatibility
  if ( !(fabs(fTR->MuD0PV[index]) < 0.02) ) return false;
  counters[MU].fill(" ... D0(pv) < 0.02");
  if ( !(fabs(fTR->MuDzPV[index]) < 1.0 ) ) return false;
  counters[MU].fill(" ... DZ(pv) < 1.0");

  // Flat isolation below 20 GeV (only for synch.: we cut at 20...)
  double hybridIso = fTR->MuRelIso03[index]*fTR->MuPt[index]/std::max((float)20.,fTR->MuPt[index]);
  if ( !(hybridIso < 0.15) ) return false;
  counters[MU].fill(" ... hybridIso < 0.15");

  if ( !(fTR->MuPtE[index]/fTR->MuPt[index] < 0.1) ) return false;
  counters[MU].fill(" ... dpt/pt < 0.1");
*/
  return true;
}


const bool JZBPFAnalysis::IsCustomPfEl(const int index){

  // kinematic acceptance
  if(!(fTR->PfEl3Pt[index]>10) )return false;
  counters[PFEL].fill(" ... PF pt > 10");
  if(!(fabs(fTR->PfEl3Eta[index]) < 2.4) ) return false;
  counters[PFEL].fill(" ... PF |eta| < 2.4");
  if(!(fTR->PfElID95[index])) return false;
  //if(!(fTR->PfElID80[index])) return false;
/*
  if ( !(fTR->ElNumberOfMissingInnerHits[index] <= 1 ) ) return false;
  counters[EL].fill(" ... missing inner hits <= 1");
  if ( !(fabs(fTR->ElD0PV[index]) < 0.04) ) return false;
  counters[EL].fill(" ... D0(pv) < 0.04");
  if ( !(fabs(fTR->ElDzPV[index]) < 1.0 ) ) return false;
  counters[EL].fill(" ... DZ(pv) < 1.0");

  // Electron ID
  int elIDWP95 = fTR->ElIDsimpleWP95relIso[index];
  if (elIDWP95!=7) return false;
  counters[EL].fill(" ... passes WP95 ID");

  // Flat isolation below 20 GeV (only for synch.)
  double hybridIso = fTR->ElRelIso03[index]
    *fTR->ElPt[index]/std::max((float)20.,fTR->ElPt[index]);
  if ( !(hybridIso < 0.15) ) return false;  
  counters[EL].fill(" ... hybridIso < 0.15");

  //   // Other choices for electron ID
  //   if ( fTR->ElIDsimpleWP90relIso[index]!=7 ) return false;
  //   counters[EL].fill("... passes WP90 ID");
  //   if ( fTR->ElIDsimpleWP80relIso[index]!=7 ) return false;
  //   counters[EL].fill("... passes WP80 ID");
  //   if ( !(fTR->ElIDMva[index]>0.4) ) return false;
  //   counters[EL].fill("... MVA>0.4");
*/
  return true;

	
}

//Get the index of the reco muon in the reco collection
int JZBPFAnalysis::getRecoMuIndex(TLorentzVector tmpPfVector) {
  int retval = -1;
  int muIndex = 0;
  float lastDeltaR=150.0;
  for(muIndex=0;muIndex<fTR->NMus;muIndex++) {
    float px= fTR->MuPx[muIndex];
    float py= fTR->MuPy[muIndex];
    float pz= fTR->MuPz[muIndex];
    float energy =  fTR->MuE[muIndex];
    TLorentzVector tmpVector(px,py,pz,energy);
    float currdeltaR = tmpPfVector.DeltaR(tmpVector);
    if(currdeltaR<lastDeltaR) {
	retval=muIndex;
	lastDeltaR=currdeltaR;
    }
  }
  return retval;
} 

//Get the index of the reco electron in the reco collection
int JZBPFAnalysis::getRecoElIndex(TLorentzVector tmpPfVector) {
  int retval = -1;
  int elIndex = 0;
  float lastDeltaR=150.0;
  for(elIndex=0;elIndex<fTR->NEles;elIndex++) {
    float px= fTR->ElPx[elIndex];
    float py= fTR->ElPy[elIndex];
    float pz= fTR->ElPz[elIndex];
    float energy =  fTR->ElE[elIndex];
    TLorentzVector tmpVector(px,py,pz,energy);
    float currdeltaR = tmpPfVector.DeltaR(tmpVector);
    if(currdeltaR<lastDeltaR) {
	retval=elIndex;
	lastDeltaR=currdeltaR;
    }
  }
  return retval;
} 





const bool JZBPFAnalysis::IsCustomMu(const int index){
  
  // Basic muon cleaning and ID

  // Acceptance cuts
  if ( !(fTR->MuPt[index] > 10) )       return false;
  counters[MU].fill(" ... pt > 10");
  if ( !(fabs(fTR->MuEta[index])<2.4) ) return false;
  counters[MU].fill(" ... |eta| < 2.4");

  // Quality cuts
  if ( !fTR->MuIsGMPT[index] )        return false;
  counters[MU].fill(" ... is global muon prompt tight");
  if ( !fTR->MuIsGlobalMuon[index] )  return false;
  counters[MU].fill(" ... is global muon");
  if ( !fTR->MuIsTrackerMuon[index] ) return false;
  counters[MU].fill(" ... is tracker muon");
  
  // Hits
  if ( !(fTR->MuNTkHits[index] >= 11) )     return false;
  counters[MU].fill(" ... nTkHits >= 11");
  if ( !(fTR->MuNPxHits[index] > 0) )       return false;
  counters[MU].fill(" ... nPxHits > 0");
  if ( !(fTR->MuNMatches[index] > 1) )      return false;
  counters[MU].fill(" ... nMatches > 1");

  // Vertex compatibility
  if ( !(fabs(fTR->MuD0PV[index]) < 0.02) ) return false;
  counters[MU].fill(" ... D0(pv) < 0.02");
  if ( !(fabs(fTR->MuDzPV[index]) < 1.0 ) ) return false;
  counters[MU].fill(" ... DZ(pv) < 1.0");

  // Flat isolation below 20 GeV (only for synch.: we cut at 20...)
  double hybridIso = fTR->MuRelIso03[index]*fTR->MuPt[index]/std::max((float)20.,fTR->MuPt[index]);
  if ( !(hybridIso < 0.15) ) return false;
  counters[MU].fill(" ... hybridIso < 0.15");

  if ( !(fTR->MuPtE[index]/fTR->MuPt[index] < 0.1) ) return false;
  counters[MU].fill(" ... dpt/pt < 0.1");

  return true;
}


const bool JZBPFAnalysis::IsCustomEl(const int index){

  // kinematic acceptance
  if(!(fTR->ElPt[index]>10) )return false;
  counters[EL].fill(" ... pt > 10");
  if(!(fabs(fTR->ElEta[index]) < 2.4) ) return false;
  counters[EL].fill(" ... |eta| < 2.4");
  if ( !(fTR->ElNumberOfMissingInnerHits[index] <= 1 ) ) return false;
  counters[EL].fill(" ... missing inner hits <= 1");
  if ( !(fabs(fTR->ElD0PV[index]) < 0.04) ) return false;
  counters[EL].fill(" ... D0(pv) < 0.04");
  if ( !(fabs(fTR->ElDzPV[index]) < 1.0 ) ) return false;
  counters[EL].fill(" ... DZ(pv) < 1.0");

  // Electron ID
  int elIDWP95 = fTR->ElIDsimpleWP95relIso[index];
  if (elIDWP95!=7) return false;
  counters[EL].fill(" ... passes WP95 ID");

  // Flat isolation below 20 GeV (only for synch.)
  double hybridIso = fTR->ElRelIso03[index]
    *fTR->ElPt[index]/std::max((float)20.,fTR->ElPt[index]);
  if ( !(hybridIso < 0.15) ) return false;  
  counters[EL].fill(" ... hybridIso < 0.15");

  //   // Other choices for electron ID
  //   if ( fTR->ElIDsimpleWP90relIso[index]!=7 ) return false;
  //   counters[EL].fill("... passes WP90 ID");
  //   if ( fTR->ElIDsimpleWP80relIso[index]!=7 ) return false;
  //   counters[EL].fill("... passes WP80 ID");
  //   if ( !(fTR->ElIDMva[index]>0.4) ) return false;
  //   counters[EL].fill("... MVA>0.4");

  return true;

	
}

const bool JZBPFAnalysis::IsCustomJet(const int index){
  // Basic Jet ID cuts (loose Jet ID)
  // See https://twiki.cern.ch/twiki/bin/view/CMS/JetID

  if ( !(fTR->CAJID_n90Hits[index] > 1) ) return false;
  counters[JE].fill(" ... n90Hits > 1");
  if ( !(fTR->CAJID_HPD[index] < 0.98)  ) return false;
  counters[JE].fill(" ... HPD < 0.98");

  if ( fabs(fTR->CAJEta[index])<2.6 ) {
    if ( !(fTR->CAJEMfrac[index] > 0.01)  ) return false;
  } else {
    if ( !(fTR->CAJEMfrac[index] > -0.9)  ) return false;
    if ( fTR->CAJPt[index] > 80 && !(fTR->CAJEMfrac[index]<1) ) return false;
  }
  counters[JE].fill(" ... pass EMfrac cut");

  return true;
}


void JZBPFAnalysis::GeneratorInfo(void) {
  // Try to find an Z->ll pair inside the acceptance
  double minPt = 20.;
  double mllCut = 20.;
  double maxEta = 2.4;
  double minJPt = 30;
  double maxJEta = 3.0;
  // First, look for leptons in acceptance
  vector<PFlepton> gLeptons;
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
        PFlepton tmpLepton;
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
  vector<PFlepton> sortedGLeptons = sortLeptonsByPt(gLeptons);

  // Store actual number of leptons passing selection
  npfEvent.genNleptons = gLeptons.size();

  // Now fill information
  TLorentzVector GenMETvector(fTR->GenMETpx,fTR->GenMETpy,0,0);
  npfEvent.genMET     = fTR->GenMET;

  // Number of good jets
  npfEvent.genNjets = 0;
  for ( int jIndex=0; jIndex<fTR->NGenJets; ++jIndex) {
    if ( fTR->GenJetPt[jIndex]<minJPt ) continue;
    if ( fabs(fTR->GenJetEta[jIndex])>maxJEta ) continue;
    ++npfEvent.genNjets;
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
      npfEvent.genPt1     = sortedGLeptons[i1].p.Pt();
      npfEvent.genId1     = sortedGLeptons[i1].type;
      npfEvent.genEta1    = sortedGLeptons[i1].p.Eta();
      npfEvent.genMID     = fTR->GenLeptonMID[sortedGLeptons[i1].index];
      npfEvent.genGMID    = fTR->GenLeptonGMID[sortedGLeptons[i1].index];
      if(sortedGLeptons.size()>1)
        {
          TLorentzVector genZvector = sortedGLeptons[i1].p + sortedGLeptons[i2].p;
          npfEvent.genRecoil  = (-GenMETvector - genZvector).Pt();
          npfEvent.genPt2     = sortedGLeptons[i2].p.Pt();
          npfEvent.genId2     = sortedGLeptons[i2].type;
          npfEvent.genEta2    = sortedGLeptons[i2].p.Eta();
          npfEvent.genZPt     = genZvector.Pt();
          npfEvent.genMll     = genZvector.M();
          npfEvent.genJZB     = npfEvent.genRecoil - genZvector.Pt();
        }
    }
}
