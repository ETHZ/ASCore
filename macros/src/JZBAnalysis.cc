#include "helper/Utilities.hh"
#include "JZBAnalysis.hh"
#include "TF1.h"
#include <time.h>
#include <TRandom.h>
#include "TF1.h"
//#include "/shome/theofil/setTDRStyle.C"

using namespace std;


#define jMax 30  // do not touch this
#define metMax 30
#define rMax 30

string sjzbversion="$Revision: 1.25 $";
string sjzbinfo="";

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


  float pt1; // leading leptons
  float pt2;

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

  int leptonPairNum;
  int leptonPairId[jMax];
  float leptonPairMass[jMax];
  float leptonPairDphi[jMax];


  int pfJetNum;
  float pfJetPt[jMax];
  float pfJetEta[jMax];
  float pfJetPhi[jMax];
  float pfJetID[jMax];
  float pfJetScale[jMax];
  float pfJetDphiMet[jMax];
  float pfHT;
  float pfGoodHT;
  float pfTightHT;

  int pfJetGoodNum;
  float pfJetGoodPt[jMax];
  float pfJetGoodEta[jMax];
  float pfJetGoodPhi[jMax];

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
  float sjzb[rMax]; // smeared JZB
  float dphi_sumJetVSZ[rMax];
  float sumJetPt[rMax];

  float weight;
  float PUweight;
  bool passed_triggers;

};

nanoEvent::nanoEvent(){};
void nanoEvent::reset()
{

  mll=0; // di-lepton system
  pt=0;
  phi=0;

  pt1=0;
  pt2=0;
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
  
  eta1=0; // leading leptons
  eta2=0;
  phi1=0;
  phi2=0;
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
  }
  pfJetGoodNum=0;
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


TTree *myTree;
TTree *myInfo;

nanoEvent nEvent;


JZBAnalysis::JZBAnalysis(TreeReader *tr, std::string dataType, bool fullCleaning) : 
  UserAnalysisBase(tr), fDataType_(dataType), fFullCleaning_(fullCleaning) {
  //	Util::SetStyle();	
  //	setTDRStyle();	
}

JZBAnalysis::~JZBAnalysis(){
}

void JZBAnalysis::Begin(){
  // Define the output file of histograms
  fHistFile = new TFile(outputFileName_.c_str(), "RECREATE");
	
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
  FILE *usernamefile;
  usernamefile = popen("whoami", "r");
  fgets(usertext, sizeof(usertext), usernamefile);
  pclose(usernamefile);
  *jzbversion=sjzbversion;
  char scmsdir[1000];
  getcwd(scmsdir,1000);
  *cmsdir=scmsdir;
  *jzbinfo=sjzbinfo;
  *user=usertext;
  time_t rawtime;
  struct tm * timeinfo;
  time (&rawtime );
  *timestamp=ctime(&rawtime);
  
  myInfo->Fill();
  myInfo->Write();

  myTree = new TTree("events","events");

  myTree->Branch("mll",&nEvent.mll,"mll/F");
  myTree->Branch("pt",&nEvent.pt,"pt/F");
  myTree->Branch("phi",&nEvent.phi,"phi/F");
  myTree->Branch("pt1",&nEvent.pt1,"pt1/F");
  myTree->Branch("pt2",&nEvent.pt2,"pt2/F");
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

  myTree->Branch("id1",&nEvent.id1,"id1/I");
  myTree->Branch("id2",&nEvent.id2,"id2/I");
  myTree->Branch("ch1",&nEvent.ch1,"ch1/I");
  myTree->Branch("ch2",&nEvent.ch2,"ch2/I");
  myTree->Branch("chid1",&nEvent.chid1,"chid1/I");
  myTree->Branch("chid2",&nEvent.chid2,"chid2/I");

  myTree->Branch("jetNum",&nEvent.jetNum,"jetNum/I");
  myTree->Branch("goodJetNum",&nEvent.goodJetNum,"goodJetNum/I");
  myTree->Branch("jetID",nEvent.jetID,"jetID[jetNum]/I");
  myTree->Branch("jetpt",nEvent.jetpt,"jetpt[jetNum]/F");
  myTree->Branch("jeteta",nEvent.jeteta,"jeteta[jetNum]/F");
  myTree->Branch("jetphi",nEvent.jetphi,"jetphi[jetNum]/F");
  myTree->Branch("jetscale",nEvent.jetscale,"jetscale[jetNum]/F");

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

  myTree->Branch("recoilpt",nEvent.recoilpt,"recoilpt[30]/F");
  myTree->Branch("dphiRecoilLep",nEvent.dphiRecoilLep,"dphiRecoilLep[30]/F");
  myTree->Branch("recoilphi",nEvent.recoilphi,"recoilphi[30]/F");
  myTree->Branch("recoileta",nEvent.recoileta,"recoileta[30]/F");
  myTree->Branch("recoilenergy",nEvent.recoilenergy,"recoilenergy[30]/F");

  myTree->Branch("vjetpt",nEvent.vjetpt,"vjetpt[30]/F");
  myTree->Branch("vjeteta",nEvent.vjeteta,"vjeteta[30]/F");
  myTree->Branch("vjetphi",nEvent.vjetphi,"vjetphi[30]/F");

  myTree->Branch("met",nEvent.met,"met[30]/F");
  myTree->Branch("metPhi",nEvent.metPhi,"metPhi[30]/F");
  myTree->Branch("dphiMetLep",nEvent.dphiMetLep,"dphiMetLep[30]/F");
  myTree->Branch("dphiMetJet",nEvent.dphiMetJet,"dphiMetJet[30]/F");
  myTree->Branch("dphiMetSumJetPt",nEvent.dphiMetSumJetPt,"dphiMetSumJetPt[30]/F");
  myTree->Branch("metPerp",nEvent.metPerp,"metPerp[30]/F");
  myTree->Branch("metPar",nEvent.metPar,"metPar[30]/F");


  myTree->Branch("eventNum",&nEvent.eventNum,"eventNum/I");
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
  myTree->Branch("pfJetID",nEvent.pfJetID,"pfJetID[pfJetNum]/F");
  myTree->Branch("pfJetScale",nEvent.pfJetScale,"pfJetScale[pfJetNum]/F");
  myTree->Branch("pfJetDphiMet",nEvent.pfJetDphiMet,"pfJetDphiMet[pfJetNum]/F");
  myTree->Branch("pfHT",&nEvent.pfHT,"pfHT/F");
  myTree->Branch("pfGoodHT",&nEvent.pfGoodHT,"pfGoodHT/F");
  myTree->Branch("pfTightHT",&nEvent.pfTightHT,"pfTightHT/F");

  myTree->Branch("pfJetGoodNum",&nEvent.pfJetGoodNum,"pfJetGoodNum/I");
  myTree->Branch("pfJetGoodPt",nEvent.pfJetGoodPt,"pfJetGoodPt[pfJetGoodNum]/F");
  myTree->Branch("pfJetGoodEta",nEvent.pfJetGoodEta,"pfJetGoodEta[pfJetGoodNum]/F");
  myTree->Branch("pfJetGoodPhi",nEvent.pfJetGoodPhi,"pfJetGoodPhi[pfJetGoodNum]/F");

  myTree->Branch("pfJetGoodNum20",&nEvent.pfJetGoodNum20,"pfJetGoodNum20/I");
  myTree->Branch("pfJetGoodNum25",&nEvent.pfJetGoodNum25,"pfJetGoodNum25/I");
  myTree->Branch("pfJetGoodNum27",&nEvent.pfJetGoodNum27,"pfJetGoodNum27/I");
  myTree->Branch("pfJetGoodNum285",&nEvent.pfJetGoodNum285,"pfJetGoodNum285/I");
  myTree->Branch("pfJetGoodNum315",&nEvent.pfJetGoodNum315,"pfJetGoodNum315/I");
  myTree->Branch("pfJetGoodNum33",&nEvent.pfJetGoodNum33,"pfJetGoodNum33/I");
  myTree->Branch("pfJetGoodNum35",&nEvent.pfJetGoodNum35,"pfJetGoodNum35/I");

  myTree->Branch("jzb",nEvent.jzb,"jzb[30]/F");
  myTree->Branch("sjzb",nEvent.sjzb,"sjzb[30]/F");
  myTree->Branch("dphi_sumJetVSZ",nEvent.dphi_sumJetVSZ,"dphi_sumJetVSZ[30]/F");
  myTree->Branch("sumJetPt",nEvent.sumJetPt,"sumJetPt[30]/F");
  myTree->Branch("weight", &nEvent.weight,"weight/F");
  myTree->Branch("PUweight",&nEvent.PUweight,"PUweight/F");
  myTree->Branch("passed_triggers", &nEvent.passed_triggers,"passed_triggers/B");


  // Define counters (so we have them in the right order)
  counters[EV].setName("Events");
  counters[TR].setName("Triggers");
  counters[MU].setName("Muons");
  counters[EL].setName("Electrons");
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
bool momentumComparator(lepton i, lepton j) { return (i.p.Pt()>j.p.Pt()); }


//------------------------------------------------------------------------------
vector<lepton> JZBAnalysis::sortLeptonsByPt(vector<lepton>& leptons) {
  
  vector<lepton> theLep = leptons;
  sort (theLep.begin(), theLep.end(), momentumComparator);
  return theLep;  
  
}


//------------------------------------------------------------------------------
const bool JZBAnalysis::passElTriggers() {
  if ( GetHLTResult("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v1") )        return true;
  if ( GetHLTResult("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v2") )        return true;
  if ( GetHLTResult("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v3") )        return true;
  return false;

}

//------------------------------------------------------------------------------
const bool JZBAnalysis::passMuTriggers() {
  if ( GetHLTResult("HLT_DoubleMu6_v1") )        return true;
  if ( GetHLTResult("HLT_DoubleMu6_v2") )        return true;
  if ( GetHLTResult("HLT_DoubleMu7_v1") )        return true;
  if ( GetHLTResult("HLT_DoubleMu7_v2") )        return true;
  return false;
} 

//______________________________________________________________________________
const bool JZBAnalysis::passEMuTriggers() {
  if ( GetHLTResult("HLT_Mu17_Ele8_CaloIdL_v1") )        return true;
  if ( GetHLTResult("HLT_Mu17_Ele8_CaloIdL_v2") )        return true;
  if ( GetHLTResult("HLT_Mu17_Ele8_CaloIdL_v3") )        return true;
  if ( GetHLTResult("HLT_Mu8_Ele17_CaloIdL_v1") )        return true;
  if ( GetHLTResult("HLT_Mu8_Ele17_CaloIdL_v2") )        return true;
  if ( GetHLTResult("HLT_Mu8_Ele17_CaloIdL_v3") )        return true;
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
      nEvent.PUweight  = GetPUWeight(fTR->PUnumInteractions);
      nEvent.weight    = GetPUWeight(fTR->PUnumInteractions);
    }

  // Trigger information
  nEvent.passed_triggers=0;
  if ( fDataType_ == "mu" || fDataType_ == "el"  || fDataType_ == "em") 
    {
      if( (fDataType_=="mu") && passMuTriggers() ) 
        {
          counters[EV].fill("... pass muon triggers");
          nEvent.passed_triggers=1;
        } 
      else if ( (fDataType_=="el") && passElTriggers() ) 
        {
          counters[EV].fill("... pass electron triggers");
          nEvent.passed_triggers=1;
        } 
      else if ( (fDataType_=="em") && passEMuTriggers() ) 
        {
          counters[EV].fill("... pass EM triggers");
          nEvent.passed_triggers=1;
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
    myTree->Fill();
    return;
  }
  counters[EV].fill("... pass good event requirements");

  vector<lepton> leptons;
  TLorentzVector genZvector; // To store the true Z vector

  // #--- muon loop
  for(int muIndex=0;muIndex<fTR->NMus;muIndex++)
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
	  
	  lepton tmpLepton;
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
	  
	  lepton tmpLepton;
	  tmpLepton.p = tmpVector;
	  tmpLepton.charge = tmpCharge;
	  tmpLepton.index = elIndex;
	  tmpLepton.type = 0;
	  tmpLepton.genPt = 0.;
	  leptons.push_back(tmpLepton);
	}
    }
  
  
  // Sort the leptons by Pt and select the two opposite-signed ones with highest Pt
  
  vector<lepton> sortedGoodLeptons = sortLeptonsByPt(leptons);

  if(sortedGoodLeptons.size() > 1) {
    
    counters[EV].fill("... has at least 2 leptons");
    int PosLepton1 = 0;
    int PosLepton2 = 1;
    
    // Check for OS combination
    for(; PosLepton2 < sortedGoodLeptons.size(); PosLepton2++) {
      if(sortedGoodLeptons[0].charge*sortedGoodLeptons[PosLepton2].charge<0) break;
    }
    if(PosLepton2 == sortedGoodLeptons.size()) {
      myTree->Fill();
      return;
    }
    counters[EV].fill("... has at least 2 OS leptons");
    
    // Preselection
    if(sortedGoodLeptons[PosLepton1].p.Pt() > 20 && sortedGoodLeptons[PosLepton2].p.Pt() > 20) {

      nEvent.eta1 = sortedGoodLeptons[PosLepton1].p.Eta();
      nEvent.pt1 = sortedGoodLeptons[PosLepton1].p.Pt();
      nEvent.phi1 = sortedGoodLeptons[PosLepton1].p.Phi();
      nEvent.ch1 = sortedGoodLeptons[PosLepton1].charge;
      nEvent.id1 = sortedGoodLeptons[PosLepton1].type; //??????
      nEvent.chid1 = (sortedGoodLeptons[PosLepton1].type+1)*sortedGoodLeptons[PosLepton1].charge;
      
      nEvent.eta2 = sortedGoodLeptons[PosLepton2].p.Eta();
      nEvent.pt2 = sortedGoodLeptons[PosLepton2].p.Pt();
      nEvent.phi2 = sortedGoodLeptons[PosLepton2].p.Phi();
      nEvent.ch2 = sortedGoodLeptons[PosLepton2].charge;
      nEvent.id2 = sortedGoodLeptons[PosLepton2].type; //??????
      nEvent.chid2 = (sortedGoodLeptons[PosLepton2].type+1)*sortedGoodLeptons[PosLepton2].charge;

      nEvent.mll=(sortedGoodLeptons[PosLepton2].p+sortedGoodLeptons[PosLepton1].p).M();
      nEvent.phi=(sortedGoodLeptons[PosLepton2].p+sortedGoodLeptons[PosLepton1].p).Phi();
      nEvent.pt=(sortedGoodLeptons[PosLepton2].p+sortedGoodLeptons[PosLepton1].p).Pt();
      nEvent.dphi=sortedGoodLeptons[PosLepton2].p.DeltaPhi(sortedGoodLeptons[PosLepton1].p);
      
    } else {
      
      //If there are less than two leptons the event is not considered
      myTree->Fill();
      return;
      
    }
    counters[EV].fill("... pass dilepton pt selection");
    
    
    // #--- construct different recoil models, initial the recoil vector will hold only the sum over the hard jets, only in the end we will add-up the lepton system
    TLorentzVector recoil(0,0,0,0); // different constructions of recoil model (under dev, need cleaning)    
    nEvent.jetNum=0;        // total jet counting
    nEvent.goodJetNum=0;    // Jets passing tighter pt cut
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

        // lepton-jet cleaning
        if ( fFullCleaning_ ) { 
          // Remove jet close to any lepton
          bool isClean(true);
          for ( size_t ilep = 0; ilep<sortedGoodLeptons.size(); ++ilep )
            if ( aJet.DeltaR(sortedGoodLeptons[ilep].p)<DRmax) isClean=false;
          if ( !isClean ) continue;
          counters[JE].fill("... pass full lepton cleaning");
        } else {
          // Remove jet close to leptons from Z candidate
          if(aJet.DeltaR(sortedGoodLeptons[PosLepton1].p)<DRmax)continue; 
          counters[JE].fill("... pass lepton 1 veto");
          if(aJet.DeltaR(sortedGoodLeptons[PosLepton2].p)<DRmax)continue;
          counters[JE].fill("... pass lepton 2 veto");
        }

        // Acceptance cuts before we use this jet
        if ( !(fabs(jeta)<2.6 && jpt>20.) ) continue;
        counters[JE].fill("... |eta|<2.6 && pt>20.");
	
        recoil+=aJet;
	
	nEvent.jetpt[nEvent.jetNum]  = aJet.Pt();
	nEvent.jeteta[nEvent.jetNum] = aJet.Eta();
	nEvent.jetphi[nEvent.jetNum] = aJet.Phi();
	nEvent.jetscale[nEvent.jetNum]  = jesC;
	if(isJetID) nEvent.jetID[nEvent.jetNum] = 1;
	
	nEvent.jetNum = nEvent.jetNum + 1 ;
	
        if ( jpt>30 ) {
          counters[JE].fill("... pt>30");
          nEvent.goodJetNum++;
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
    nEvent.pfJetNum=0;
    nEvent.pfJetGoodNum=0;
    nEvent.pfJetGoodNum20=0;
    nEvent.pfJetGoodNum25=0;
    nEvent.pfJetGoodNum27=0;
    nEvent.pfJetGoodNum285=0;
    nEvent.pfJetGoodNum315=0;
    nEvent.pfJetGoodNum33=0;
    nEvent.pfJetGoodNum35=0;

    vector<lepton> pfGoodJets;
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
          if(aJet.DeltaR(sortedGoodLeptons[PosLepton1].p)<DRmax)continue; 
          counters[PJ].fill("... pass lepton 1 veto");
          if(aJet.DeltaR(sortedGoodLeptons[PosLepton2].p)<DRmax)continue;
          counters[PJ].fill("... pass lepton 2 veto");
        }
	

        // Keep jets over min. pt threshold
        if ( !(jpt>20) ) continue;  
        counters[PJ].fill("... pt>20.");

	nEvent.pfJetPt[nEvent.pfJetNum] = jpt;
	nEvent.pfJetEta[nEvent.pfJetNum] = jeta;
	nEvent.pfJetPhi[nEvent.pfJetNum] = jphi;
        nEvent.pfJetScale[nEvent.pfJetNum] = jesC;
        nEvent.pfJetID[nEvent.pfJetNum]  = static_cast<int>(isJetID);
	nEvent.pfJetDphiMet[nEvent.pfJetNum] = aJet.DeltaPhi(pfMETvector);
	nEvent.pfJetNum = nEvent.pfJetNum +1;
	nEvent.pfHT += jpt;

        // Keep central jets
        if ( !(fabs(jeta)<2.6 ) ) continue;
        counters[PJ].fill("... |eta|<2.6");

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
          nEvent.pfJetGoodNum++;
        }
        if ( jpt>20. ) nEvent.pfJetGoodNum20++;
        if ( jpt>25. ) nEvent.pfJetGoodNum25++;
        if ( jpt>27. ) nEvent.pfJetGoodNum27++;
        if ( jpt>28.5 ) nEvent.pfJetGoodNum285++;
        if ( jpt>31.5 ) nEvent.pfJetGoodNum315++;
        if ( jpt>33. ) nEvent.pfJetGoodNum33++;
        if ( jpt>35. ) nEvent.pfJetGoodNum35++;
      }
    
    
    int index;
    if(recoil.Pt()!=0) // so far we had not added the lepton system in the recoil, so our recoil represents the sumJPt (ugly but it should work)
      {
	index=0;
	nEvent.vjetpt[index]=recoil.Pt();  // vjet = vector sum of jets, vjetpt = sumJPt
	nEvent.vjeteta[index]=recoil.Eta();
	nEvent.vjetphi[index]=recoil.Phi();
      }
    
    
    TLorentzVector s1 = sortedGoodLeptons[PosLepton1].p;
    TLorentzVector s2 = sortedGoodLeptons[PosLepton2].p;

    nEvent.met[0]=fTR->RawMET;
    nEvent.met[1]=0.; // Not there anymore: fTR->MuJESCorrMET;
    nEvent.met[2]=fTR->TCMET;
    nEvent.met[3]=fTR->MuJESCorrMET;
    nEvent.met[4]=fTR->PFMET;
    nEvent.met[5]=fTR->SumEt;

    TLorentzVector caloVector(0,0,0,0); // for constructing SumJPt from raw calomet
    TLorentzVector pfJetVector(0,0,0,0); // for constructing SumJPt from pf jets, as Pablo
    TLorentzVector pfNoCutsJetVector(0,0,0,0); // for constructing SumJPt from pfmet (unclustered), as Kostas
    TLorentzVector tcNoCutsJetVector(0,0,0,0); // for constructing SumJPt from tcmet (unclustered), new
    nEvent.metPhi[0]=caloMETvector.Phi();
    nEvent.metPhi[1]=0.;
    nEvent.metPhi[2]=tcMETvector.Phi();
    nEvent.metPhi[3]=0.;
    nEvent.metPhi[4]=pfMETvector.Phi();
    nEvent.metPhi[5]=0.;
    
    
    // Remove electrons from MET
    caloVector = -caloMETvector;
    if ( sortedGoodLeptons[PosLepton1].type == 0 ) caloVector -= s1;
    if ( sortedGoodLeptons[PosLepton2].type == 0 ) caloVector -= s2;

    // remove the leptons from PFMET and tcMET blublu
    pfNoCutsJetVector = -pfMETvector - s1 - s2;
    tcNoCutsJetVector = -tcMETvector - s1 - s2;

    // #--- different versions of JZB
    nEvent.dphi_sumJetVSZ[0]=caloVector.DeltaPhi(s1+s2); // DPhi between Z and SumJpt
    nEvent.sumJetPt[0]=caloVector.Pt();
    nEvent.jzb[0] = caloVector.Pt() - (s1+s2).Pt(); // calib issue of rawcalomet wrt lepton energy scale, under develop
    
    nEvent.dphi_sumJetVSZ[1] = pfNoCutsJetVector.DeltaPhi(s1+s2); 
    nEvent.sumJetPt[1] = pfNoCutsJetVector.Pt(); 
    nEvent.jzb[1] = pfNoCutsJetVector.Pt() - (s1+s2).Pt(); // to be used with pfMET
    nEvent.sjzb[1] = GausRandom(nEvent.jzb[1]+1.3,7); // to be used with pfMET

    nEvent.dphi_sumJetVSZ[2] = recoil.DeltaPhi(s1+s2);  // recoil is not yet a recoil but the sumJPt, since the leptons will be added only later (ugly)
    nEvent.sumJetPt[2] = recoil.Pt(); 
    nEvent.jzb[2] = recoil.Pt() - (s1+s2).Pt(); // to be used recoil met (recoilpt[0])    
    nEvent.jzb[3] = sumOfPFJets.Pt() - (s1+s2).Pt(); // to be used recoil met (recoilpt[0])
    nEvent.sumJetPt[3] = sumOfPFJets.Pt();

    nEvent.dphi_sumJetVSZ[4] = tcNoCutsJetVector.DeltaPhi(s1+s2); // tcJZB
    nEvent.sumJetPt[4] = tcNoCutsJetVector.Pt(); 
    nEvent.jzb[4] = tcNoCutsJetVector.Pt() - (s1+s2).Pt(); // to be used with tcMET

    // --- recoil met and pf recoil met
    nEvent.met[6] = (sumOfPFJets + s1 + s2).Pt(); 
    nEvent.met[7] = (recoil + s1 + s2).Pt();
    
    // ----------------------------------------
    recoil+=s1+s2;   // now add also the leptons to the recoil! to form the complete recoil model

    if(recoil.Pt()!=0)
      {
	index=0;
	nEvent.recoilpt[index]=recoil.Pt();
	nEvent.recoileta[index]=recoil.Eta();
	nEvent.recoilphi[index]=recoil.Phi();
      }
    
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
    nEvent.dphitcMet1 = sortedGoodLeptons[PosLepton1].p.DeltaPhi(tcMETvector);
    nEvent.dphitcMet2 = sortedGoodLeptons[PosLepton2].p.DeltaPhi(tcMETvector);
    nEvent.dphipfRecoilMet1 = sortedGoodLeptons[PosLepton1].p.DeltaPhi(-sumOfPFJets - s1 - s2); // pf recoil met
    nEvent.dphipfRecoilMet2 = sortedGoodLeptons[PosLepton2].p.DeltaPhi(-sumOfPFJets - s1 - s2); // pf recoil met
    
    // Store minimum dphi between some mets and any kind of lepton
    for ( size_t i=0; i<sortedGoodLeptons.size(); ++i ) {
      TLorentzVector lp(sortedGoodLeptons[i].p);
      if ( fabs(recoil.DeltaPhi(lp))>nEvent.dphiRecoilLep[0] )
        nEvent.dphiRecoilLep[0] = recoil.DeltaPhi(lp);
      if ( fabs(pfMETvector.DeltaPhi(lp))>nEvent.dphiMetLep[4] )
        nEvent.dphiMetLep[4] = pfMETvector.DeltaPhi(lp);
      if ( fabs((sumOfPFJets + s1 + s2).DeltaPhi(lp))> nEvent.dphiMetLep[6] )
        nEvent.dphiMetLep[6] = (sumOfPFJets + s1 + s2).DeltaPhi(lp);
      if ( fabs((recoil + s1 + s2).DeltaPhi(lp)) > nEvent.dphiMetLep[7] )
        nEvent.dphiMetLep[7] = (recoil + s1 + s2).DeltaPhi(lp);
    }

    // Store minimum dphi between some mets and any good jet
    for ( size_t i=0; i<pfGoodJets.size(); ++i ) {
      TLorentzVector jp(pfGoodJets[i].p);
      if ( fabs(pfMETvector.DeltaPhi(jp))>nEvent.dphiMetJet[4] )
        nEvent.dphiMetJet[4] = pfMETvector.DeltaPhi(jp);
    }
    nEvent.dphiMetSumJetPt[4] = pfNoCutsJetVector.DeltaPhi(pfMETvector);

    // Store some additional MET information
    nEvent.metPar[4]  = pfMETvector.Dot(s1+s2);
    nEvent.metPerp[4] = pfMETvector.Perp((s1+s2).Vect());
    
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

      nEvent.genMID1     = fTR->GenLeptonMID[i1]; // WW study
      nEvent.genMID2     = fTR->GenLeptonMID[i2]; // WW study

      if(sortedGoodLeptons.size()>=3) {i3=sortedGoodLeptons[2].index;nEvent.genMID3 = fTR->GenLeptonMID[i3];} else nEvent.genMID3=-999; // WW study
      if(sortedGoodLeptons.size()>=4) {i4=sortedGoodLeptons[3].index;nEvent.genMID4 = fTR->GenLeptonMID[i4];} else nEvent.genMID4=-999; // WW study
      if(sortedGoodLeptons.size()>=5) {i5=sortedGoodLeptons[4].index;nEvent.genMID5 = fTR->GenLeptonMID[i5];} else nEvent.genMID5=-999; // WW study

      nEvent.genGMID1    = fTR->GenLeptonGMID[i1]; // WW study
      nEvent.genGMID2    = fTR->GenLeptonGMID[i2]; // WW study

      nEvent.genPt1Sel    = genLep1.Pt();
      nEvent.genPt2Sel    = genLep2.Pt();
      nEvent.genEta1Sel   = genLep1.Eta();
      nEvent.genEta2Sel   = genLep2.Eta();
      nEvent.genId1Sel    = (sortedGoodLeptons[PosLepton1].type?fTR->MuGenID[i1]:fTR->ElGenID[i1]);
      nEvent.genId2Sel    = (sortedGoodLeptons[PosLepton2].type?fTR->MuGenID[i2]:fTR->ElGenID[i2]);
      nEvent.genJZBSel    = nEvent.genRecoilSel - (genLep1 + genLep2).Pt();
	


      TLorentzVector tmpVector;
      tmpVector.SetPtEtaPhiM(fTR->GenLeptonPt[i1],fTR->GenLeptonEta[i1],fTR->GenLeptonPhi[i1],0.);
      TLorentzVector l1=tmpVector;
      tmpVector.SetPtEtaPhiM(fTR->GenLeptonPt[i2],fTR->GenLeptonEta[i2],fTR->GenLeptonPhi[i2],0.);
      TLorentzVector l2=tmpVector;
      TLorentzVector genZvector = l1 + l2;
      nEvent.genRecoil  = (-GenMETvector - genZvector).Pt();
      nEvent.genZPt     = genZvector.Pt();
      nEvent.genMll     = genZvector.M();
      nEvent.genJZB     = nEvent.genRecoil - genZvector.Pt();
      nEvent.genPt1     = l1.Pt();
      nEvent.genPt2     = l2.Pt();
      nEvent.genId1     = fTR->GenLeptonID[i1];
      nEvent.genId2     = fTR->GenLeptonID[i2];
      nEvent.genEta1    = fTR->GenLeptonEta[i1];
      nEvent.genEta2    = fTR->GenLeptonEta[i2];
    }
    myTree->Fill();
  }
    
}

void JZBAnalysis::End(){
  fHistFile->cd();	

  myTree->Write();

  // Dump statistics
  if (1) { // Put that to 0 if you are annoyed
    std::cout << setfill('=') << std::setw(70) << "" << std::endl;
    std::cout << "Statistics" << std::endl;
    std::cout << setfill('-') << std::setw(70) << "" << setfill(' ') << std::endl;
    for ( counters_t iCount=count_begin; iCount<count_end; 
          iCount = counters_t(iCount+1) ) {
      counters[iCount].print();
      std::cout << setfill('-') << std::setw(70) << "" << setfill(' ') << std::endl;    
    }
  }

  fHistFile->Close();
}

template<class T>
std::string JZBAnalysis::any2string(T i)
{
  std::ostringstream buffer;
  buffer << i;
  return buffer.str();
}


const bool JZBAnalysis::IsCustomMu(const int index){
  
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


const bool JZBAnalysis::IsCustomEl(const int index){

  // kinematic acceptance
  if(!(fTR->ElPt[index]>10) )return false;
  counters[EL].fill(" ... pt > 10");
  if(!(fabs(fTR->ElEta[index]) < 2.5) ) return false;
  counters[EL].fill(" ... |eta| < 2.5");
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

const bool JZBAnalysis::IsCustomJet(const int index){
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


void JZBAnalysis::GeneratorInfo(void) {
  // Try to find an Z->ll pair inside the acceptance
  double minPt = 20.;
  double mllCut = 20.;
  double maxEtaEl = 2.5, maxEtaMu = 2.4;
  // First, look for leptons in acceptance and coming from Z
  vector<lepton> gLeptons;
  for ( int gIndex=0;gIndex<fTR->NGenLeptons; ++gIndex ) {
    if ( fTR->GenLeptonPt[gIndex]>minPt && 
         (( abs(fTR->GenLeptonID[gIndex])==11 
            && fabs(fTR->GenLeptonEta[gIndex])<maxEtaEl) 
          || ( abs(fTR->GenLeptonID[gIndex])==13 
               && fabs(fTR->GenLeptonEta[gIndex])<maxEtaMu) ) )
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
        //if ( fTR->GenLeptonMID[gIndex] ==23 ) gLeptons.push_back(tmpLepton); // WW study
      }         
  }
  // Stop here if not even two leptons from Z
  //  if ( gLeptons.size()<2 ) return;
  // Problem: we have too many leptons from Z... Will take two first
  if ( gLeptons.size()>2 ) 
    std::cerr << "*** WARNING: more than one Z found: " 
              << gLeptons.size() << std::endl;

  // Check for Z in the mll window

  //if ( fabs(genZvector.M()-91.2)>mllCut ) return; // WW study


  // Now fill information
  TLorentzVector GenMETvector(fTR->GenMETpx,fTR->GenMETpy,0,0);
  
  nEvent.genMET     = fTR->GenMET;
  
  if(gLeptons.size()>0)
    {
      nEvent.genPt1     = gLeptons[0].p.Pt();
      nEvent.genId1     = gLeptons[0].type;
      nEvent.genEta1    = gLeptons[0].p.Eta();
      nEvent.genMID     = fTR->GenLeptonMID[gLeptons[0].index];
      nEvent.genGMID    = fTR->GenLeptonGMID[gLeptons[0].index];
      if(gLeptons.size()>1)
        {
          TLorentzVector genZvector = gLeptons[0].p + gLeptons[1].p;
          nEvent.genRecoil  = (-GenMETvector - genZvector).Pt();
          nEvent.genPt2     = gLeptons[1].p.Pt();
          nEvent.genId2     = gLeptons[1].type;
          nEvent.genEta2    = gLeptons[1].p.Eta();
          nEvent.genZPt     = genZvector.Pt();
          nEvent.genMll     = genZvector.M();
          nEvent.genJZB     = nEvent.genRecoil - genZvector.Pt();
        }
    }
}
