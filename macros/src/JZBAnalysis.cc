#include "helper/Utilities.hh"
#include "JZBAnalysis.hh"
//#include "/shome/theofil/setTDRStyle.C"

using namespace std;


#define jMax 30  // do not touch this
#define metMax 30
#define rMax 30



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
  float eta1; // leading leptons
  float eta2;
  float phi1;
  float phi2;
  float dphi;
  int id1;
  int id2;
  int iso1;
  int iso2;
  int ch1;
  int ch2;


  int centraljetNum;
  float centraljetpt[jMax]; // jets in the barrel
  float centraljeteta[jMax];
  float centraljetphi[jMax];
  float centraljetpx[jMax];
  float centraljetpy[jMax];
  float centraljetpz[jMax];
  float centraljetenergy[jMax];
  float centraljetscale[jMax];
  int centraljetID[jMax];

  int jetNum;
  int goodJetNum;
  float jetpt[jMax]; // jets in barrel + endcaps
  float jeteta[jMax];
  float jetphi[jMax];
  float jetpx[jMax];
  float jetpy[jMax];
  float jetpz[jMax];
  float jetenergy[jMax];
  float jetscale[jMax];
  int jetID[jMax];


  int pfJetNum;
  float pfJetPt[jMax];
  float pfJetEta[jMax];
  float pfJetPhi[jMax];
  float pfHT;
  float pfGoodHT;
  float pfTightHT;

  int pfJetGoodNum;
  float pfJetGoodPt[jMax];
  float pfJetGoodEta[jMax];
  float pfJetGoodPhi[jMax];

  float recoilpt[rMax];
  float dphiRecoilLep[rMax];
  float vjetpt[rMax];
  float vjeteta[rMax];
  float vjetphi[rMax];
  float recoilenergy[rMax];
  float recoilphi[rMax];
  float recoileta[rMax];

  float met[metMax];
  float dphiMetLep[metMax];
  int eventNum;
  int runNum;
  int lumi;
  int goodVtx;
  int numVtx;
  float totEvents; // tot events processed by the ntuple producer (job submission efficiency), no need to keep this as int, better as float

  float jzb[rMax];
  float dphi_sumJetVSZ[rMax];
  float sumJetPt[rMax];
};

nanoEvent::nanoEvent(){};
void nanoEvent::reset()
{

  mll=0; // di-lepton system
  pt=0;
  phi=0;

  eta1=0; // leading leptons
  eta2=0;
  phi1=0;
  phi2=0;
  dphi=0;
  id1=0;
  id2=0;
  iso1=0;
  iso2=0;
  ch1=0;
  ch2=0;


  for(int jCounter=0;jCounter<jMax;jCounter++){
    centraljetpt[jCounter]=0; // jets in the barrel
    centraljeteta[jCounter]=0;
    centraljetphi[jCounter]=0;
    centraljetpx[jCounter]=0;
    centraljetpy[jCounter]=0;
    centraljetpz[jCounter]=0;
    centraljetenergy[jCounter]=0;
    centraljetscale[jCounter]=0;
    centraljetID[jCounter]=0;
  }
  centraljetNum=0;

  for(int jCounter=0;jCounter<jMax;jCounter++){
    jetpt[jCounter]=0; // jets in barrel + endcaps
    jeteta[jCounter]=0;
    jetphi[jCounter]=0;
    jetpx[jCounter]=0;
    jetpy[jCounter]=0;
    jetpz[jCounter]=0;
    jetenergy[jCounter]=0;
    jetscale[jCounter]=0;
    jetID[jCounter]=0;
  }
  jetNum=0;
  goodJetNum=0;

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
    dphiMetLep[metCounter]=0;
  }

  for(int jCounter=0;jCounter<jMax;jCounter++){
    pfJetPt[jCounter]=0;
    pfJetEta[jCounter]=0;
    pfJetPhi[jCounter]=0;
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

  eventNum=0;
  runNum=0;
  lumi=0;
  goodVtx=0;
  numVtx=0;
  totEvents=0;

  for(int rCounter=0;rCounter<rMax;rCounter++){
    jzb[rCounter]=0;
    dphi_sumJetVSZ[rCounter]=0;
    sumJetPt[rCounter]=0;
  }
}


TTree *myTree;
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

  myTree = new TTree("events","events");

  myTree->Branch("mll",&nEvent.mll,"mll/F");
  myTree->Branch("pt",&nEvent.pt,"pt/F");
  myTree->Branch("phi",&nEvent.phi,"phi/F");
  myTree->Branch("pt1",&nEvent.pt1,"pt1/F");
  myTree->Branch("pt2",&nEvent.pt2,"pt2/F");
  myTree->Branch("eta1",&nEvent.eta1,"eta1/F");
  myTree->Branch("eta2",&nEvent.eta2,"eta2/F");
  myTree->Branch("phi1",&nEvent.phi1,"phi1/F");
  myTree->Branch("phi2",&nEvent.phi2,"phi2/F");
  myTree->Branch("dphi",&nEvent.dphi,"dphi/F");

  myTree->Branch("id1",&nEvent.id1,"id1/I");
  myTree->Branch("id2",&nEvent.id2,"id2/I");
  myTree->Branch("iso1",&nEvent.iso1,"iso1/I");
  myTree->Branch("iso2",&nEvent.iso2,"iso2/I");
  myTree->Branch("ch1",&nEvent.ch1,"ch1/I");
  myTree->Branch("ch2",&nEvent.ch2,"ch2/I");

  myTree->Branch("centraljetNum",&nEvent.centraljetNum,"centraljetNum/I");
  myTree->Branch("centraljetID",nEvent.centraljetID,"centraljetID[centraljetNum]/I");
  myTree->Branch("centraljetpt",nEvent.centraljetpt,"centraljetpt[centraljetNum]/F");
  myTree->Branch("centraljeteta",nEvent.centraljeteta,"centraljeteta[centraljetNum]/F");
  myTree->Branch("centraljetphi",nEvent.centraljetphi,"centraljetphi[centraljetNum]/F");
  myTree->Branch("centraljetpx",nEvent.centraljetpx,"centraljetpx[centraljetNum]/F");
  myTree->Branch("centraljetpy",nEvent.centraljetpy,"centraljetpy[centraljetNum]/F");
  myTree->Branch("centraljetpz",nEvent.centraljetpz,"centraljetpz[centraljetNum]/F");
  myTree->Branch("centraljetenergy",nEvent.centraljetenergy,"centraljetenergy[centraljetNum]/F");
  myTree->Branch("centraljetscale",nEvent.centraljetscale,"centraljetscale[centraljetNum]/F");

  myTree->Branch("jetNum",&nEvent.jetNum,"jetNum/I");
  myTree->Branch("goodJetNum",&nEvent.goodJetNum,"goodJetNum/I");
  myTree->Branch("jetID",nEvent.jetID,"jetID[jetNum]/I");
  myTree->Branch("jetpt",nEvent.jetpt,"jetpt[jetNum]/F");
  myTree->Branch("jeteta",nEvent.jeteta,"jeteta[jetNum]/F");
  myTree->Branch("jetphi",nEvent.jetphi,"jetphi[jetNum]/F");
  myTree->Branch("jetpx",nEvent.jetpx,"jetpx[jetNum]/F");
  myTree->Branch("jetpy",nEvent.jetpy,"jetpy[jetNum]/F");
  myTree->Branch("jetpz",nEvent.jetpz,"jetpz[jetNum]/F");
  myTree->Branch("jetenergy",nEvent.jetenergy,"jetenergy[jetNum]/F");
  myTree->Branch("jetscale",nEvent.jetscale,"jetscale[jetNum]/F");

  myTree->Branch("recoilpt",nEvent.recoilpt,"recoilpt[30]/F");
  myTree->Branch("dphiRecoilLep",nEvent.dphiRecoilLep,"dphiRecoilLep[30]/F");
  myTree->Branch("recoilphi",nEvent.recoilphi,"recoilphi[30]/F");
  myTree->Branch("recoileta",nEvent.recoileta,"recoileta[30]/F");
  myTree->Branch("recoilenergy",nEvent.recoilenergy,"recoilenergy[30]/F");

  myTree->Branch("vjetpt",nEvent.vjetpt,"vjetpt[30]/F");
  myTree->Branch("vjeteta",nEvent.vjeteta,"vjeteta[30]/F");
  myTree->Branch("vjetphi",nEvent.vjetphi,"vjetphi[30]/F");

  myTree->Branch("met",nEvent.met,"met[30]/F");
  myTree->Branch("dphiMetLep",nEvent.dphiMetLep,"dphiMetLep[30]/F");


  myTree->Branch("eventNum",&nEvent.eventNum,"eventNum/I");
  myTree->Branch("runNum",&nEvent.runNum,"runNum/I");
  myTree->Branch("lumi",&nEvent.lumi,"lumi/I");
  myTree->Branch("goodVtx",&nEvent.goodVtx,"goodVtx/I");
  myTree->Branch("numVtx",&nEvent.numVtx,"numVtx/I");
  myTree->Branch("totEvents",&nEvent.totEvents,"totEvents/F");

  myTree->Branch("pfJetNum",&nEvent.pfJetNum,"pfJetNum/I");
  myTree->Branch("pfJetPt",nEvent.pfJetPt,"pfJetPt[pfJetNum]/F");
  myTree->Branch("pfJetEta",nEvent.pfJetEta,"pfJetEta[pfJetNum]/F");
  myTree->Branch("pfJetPhi",nEvent.pfJetPhi,"pfJetPhi[pfJetNum]/F");
  myTree->Branch("pfHT",&nEvent.pfHT,"pfHT/F");
  myTree->Branch("pfGoodHT",&nEvent.pfGoodHT,"pfGoodHT/F");
  myTree->Branch("pfTightHT",&nEvent.pfTightHT,"pfTightHT/F");

  myTree->Branch("pfJetGoodNum",&nEvent.pfJetGoodNum,"pfJetGoodNum/I");
  myTree->Branch("pfJetGoodPt",nEvent.pfJetGoodPt,"pfJetGoodPt[pfJetGoodNum]/F");
  myTree->Branch("pfJetGoodEta",nEvent.pfJetGoodEta,"pfJetGoodEta[pfJetGoodNum]/F");
  myTree->Branch("pfJetGoodPhi",nEvent.pfJetGoodPhi,"pfJetGoodPhi[pfJetGoodNum]/F");

  myTree->Branch("jzb",nEvent.jzb,"jzb[30]/F");
  myTree->Branch("dphi_sumJetVSZ",nEvent.dphi_sumJetVSZ,"dphi_sumJetVSZ[30]/F");
  myTree->Branch("sumJetPt",nEvent.sumJetPt,"sumJetPt[30]/F");

  // Define event counters (so we have them in the right order)
  counters[EV].fill("All events",0.);
  counters[EV].fill("... pass electron triggers",0.);
  counters[EV].fill("... pass muon triggers",0.);
  counters[EV].fill("... pass all trigger requirements",0.);
  std::string types[3] = { "ee","mm","em" };
  for ( size_t itype=0; itype<3; ++itype ) {
    counters[EV].fill("... "+types[itype]+" pairs",0.); 
    counters[EV].fill("... "+types[itype]+" + 2 jets",0.);
    counters[EV].fill("... "+types[itype]+" + 2 jets + PFMET>50 + HT>100 + Z veto",0.);
    counters[EV].fill("... "+types[itype]+" + 2 jets + PFMET>50 + HT>100 + require Z",0.);
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

  if ( GetHLTResult("HLT_Ele17_SW_TightEleId_L1R") )                return true;
  if ( GetHLTResult("HLT_Ele17_SW_TighterEleId_L1R_v1") )           return true;
  if ( GetHLTResult("HLT_DoubleEle15_SW_L1R_v1") )                  return true;
  if ( GetHLTResult("HLT_DoubleEle17_SW_L1R_v1") )                  return true;
  if ( GetHLTResult("HLT_Ele17_SW_TightCaloEleId_Ele8HE_L1R_v1") )  return true;
  if ( GetHLTResult("HLT_Ele17_SW_TightCaloEleId_SC8HE_L1R_v1") )   return true;
  if ( GetHLTResult("HLT_DoubleEle10_SW_L1R") )                     return true;
  if ( GetHLTResult("HLT_DoubleEle5_SW_L1R") )                      return true;
  if ( GetHLTResult("HLT_Ele17_SW_CaloEleId_L1R") )                 return true;
  if ( GetHLTResult("HLT_Ele17_SW_EleId_L1R") )                     return true;
  if ( GetHLTResult("HLT_Ele17_SW_LooseEleId_L1R") )                return true;
  if ( GetHLTResult("HLT_Ele15_SW_CaloEleId_L1R") )                 return true;
  if ( GetHLTResult("HLT_Ele15_SW_EleId_L1R") )                     return true;
  if ( GetHLTResult("HLT_Ele15_SW_L1R") )                           return true;
  if ( GetHLTResult("HLT_Ele15_LW_L1R") )                           return true;
  if ( GetHLTResult("HLT_Ele20_SW_L1R") )                           return true;
  if ( GetHLTResult("HLT_Ele10_SW_EleId_L1R") )                     return true;
  if ( GetHLTResult("HLT_Ele10_LW_EleId_L1R") )                     return true;
  if ( GetHLTResult("HLT_Ele10_LW_L1R") )                           return true;
  if ( GetHLTResult("HLT_Ele10_SW_L1R") )                           return true;
  if ( GetHLTResult("HLT_Ele17_SW_TighterEleIdIsol_L1R_v2") )       return true;
  if ( GetHLTResult("HLT_Ele22_SW_TighterEleId_L1R_v2") )           return true;
  if ( GetHLTResult("HLT_Ele32_SW_TightCaloEleIdTrack_L1R_v1") )    return true;
  if ( GetHLTResult("HLT_Ele32_SW_TighterEleId_L1R_v2") )           return true;
  if ( GetHLTResult("HLT_Ele27_SW_TightCaloEleIdTrack_L1R_v1") )    return true;
  if ( GetHLTResult("HLT_Ele22_SW_TighterCaloIdIsol_L1R_v2") )      return true;
  if ( GetHLTResult("HLT_Ele22_SW_TighterEleId_L1R_v3") )           return true;
  if ( GetHLTResult("HLT_Ele22_SW_TighterCaloIdIsol_L1R_v2") )      return true;

  return false;

}

//------------------------------------------------------------------------------
const bool JZBAnalysis::passMuTriggers() {

  if ( GetHLTResult("HLT_Mu9") )            return true;
  if ( GetHLTResult("HLT_Mu11") )           return true;
  if ( GetHLTResult("HLT_Mu15") )           return true;
  if ( GetHLTResult("HLT_Mu9_v1") )         return true;
  if ( GetHLTResult("HLT_Mu11_v1") )        return true;
  if ( GetHLTResult("HLT_Mu15_v1") )        return true;
  if ( GetHLTResult("HLT_DoubleMu3") )      return true;
  if ( GetHLTResult("HLT_DoubleMu3_v2") )   return true;
  if ( GetHLTResult("HLT_DoubleMu5_v1") )   return true;
  if ( GetHLTResult("HLT_Mu5_Ele5_v1") )    return true;
  if ( GetHLTResult("HLT_Mu5_Ele9_v1") )    return true;
  if ( GetHLTResult("HLT_Mu11_Ele8_v1") )   return true;
  if ( GetHLTResult("HLT_Mu8_Ele8_v1") )    return true;
  if ( GetHLTResult("HLT_Mu5_Ele13_v2") )   return true;
  if ( GetHLTResult("HLT_Mu5_Ele13_v2") )   return true;
  if ( GetHLTResult("HLT_Mu5_Ele17_v1") )   return true;
  if ( GetHLTResult("HLT_Mu7") )            return true;
  if ( GetHLTResult("HLT_Mu5") )            return true;
  if ( GetHLTResult("HLT_Mu17_v1") )        return true;
  if ( GetHLTResult("HLT_Mu19_v1") )        return true; 
  return false;
    
} 




void JZBAnalysis::Analyze(){
  
  counters[EV].fill("All events");
  if ( fDataType_ != "mc" ) {
    //FIXME: NEED TO TEST IN WHICH STREAM WE ARE (electron/muon)
    if( (fDataType_=="mu") && passMuTriggers() ) {
      counters[EV].fill("... pass muon triggers");
    } else if ( (fDataType_=="el") && passElTriggers() ) {
      counters[EV].fill("... pass electron triggers");
    } else {
      return;
    }
    counters[EV].fill("... pass all trigger requirements");
  }
  
  // #--- analysis global parameters
  double DRmax=0.4; // veto jets in a cone of DRmax close to the lepton
  nEvent.reset();
  
  // #--- Vertex info
  nEvent.numVtx = fTR->NVrtx;
  float rho = sqrt(fTR->PrimVtxx*fTR->PrimVtxx + fTR->PrimVtxy*fTR->PrimVtxy);
  if(fTR->PrimVtxGood) nEvent.goodVtx |=2; // save bits of vertex quality
  if (fTR->PrimVtxGood==0 && fTR->PrimVtxIsFake==0 && fTR->PrimVtxNdof>4 && fabs(fTR->PrimVtxz)<24 && rho<2)
    nEvent.goodVtx |=4;
  
  // Good event requirement: essentially vertex requirements
  if ( !IsGoodEvent() ) return;
  counters[EV].fill("... pass good event requirements");

  
  vector<lepton> leptons;
  
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

	  leptons.push_back(tmpLepton);
	}
    }
  
  
  // Sort the leptons by Pt and select the two opposite-signed ones with highest Pt
  
  nEvent.reset();
  
  vector<lepton> sortedGoodLeptons = sortLeptonsByPt(leptons);
  
  if(sortedGoodLeptons.size() > 1) {
    
    counters[EV].fill("... has at least 2 leptons");
    int PosLepton1 = 0;
    int PosLepton2 = 1;
    
    // Check for OS combination
    for(; PosLepton2 < sortedGoodLeptons.size(); PosLepton2++) {
      if(sortedGoodLeptons[0].charge*sortedGoodLeptons[PosLepton2].charge<0) break;
    }
    if(PosLepton2 == sortedGoodLeptons.size())return;
    counters[EV].fill("... has at least 2 OS leptons");
    
    // Preselection
    if(sortedGoodLeptons[PosLepton1].p.Pt() > 20 && sortedGoodLeptons[PosLepton2].p.Pt() > 20) {
      
      nEvent.eta1 = sortedGoodLeptons[PosLepton1].p.Eta();
      nEvent.pt1 = sortedGoodLeptons[PosLepton1].p.Pt();
      nEvent.phi1 = sortedGoodLeptons[PosLepton1].p.Phi();
      nEvent.ch1 = sortedGoodLeptons[PosLepton1].charge;
      nEvent.id1 = sortedGoodLeptons[PosLepton1].type; //??????
      
      nEvent.eta2 = sortedGoodLeptons[PosLepton2].p.Eta();
      nEvent.pt2 = sortedGoodLeptons[PosLepton2].p.Pt();
      nEvent.phi2 = sortedGoodLeptons[PosLepton2].p.Phi();
      nEvent.ch2 = sortedGoodLeptons[PosLepton2].charge;
      nEvent.id2 = sortedGoodLeptons[PosLepton2].type; //??????
      
      nEvent.mll=(sortedGoodLeptons[PosLepton2].p+sortedGoodLeptons[PosLepton1].p).M();
      nEvent.phi=(sortedGoodLeptons[PosLepton2].p+sortedGoodLeptons[PosLepton1].p).Phi();
      nEvent.pt=(sortedGoodLeptons[PosLepton2].p+sortedGoodLeptons[PosLepton1].p).Pt();
      nEvent.dphi=sortedGoodLeptons[PosLepton2].p.DeltaPhi(sortedGoodLeptons[PosLepton2].p);
      
    } else {
      
      //If there are less than two leptons the event is not considered
      return;
      
    }
    counters[EV].fill("... pass dilepton pt selection");
    
    
    // #--- construct different recoil models, initial the recoil vector will hold only the sum over the hard jets, only in the end we will add-up the lepton system
    TLorentzVector recoil(0,0,0,0); // different constructions of recoil model (under dev, need cleaning)    
    nEvent.centraljetNum=0; // barrel jet counting
    nEvent.jetNum=0;        // total jet counting
    nEvent.goodJetNum=0;    // Jets passing tighter pt cut
    for(int i =0 ; i<fTR->NJets;i++) // jet loop
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
	
	if (!isJetID) continue; //consider only Jets passing JetID
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
        if ( fabs(jeta)>2.6 && jpt>15. ) continue;
        counters[JE].fill("... |eta|<2.6 && pt>15.");
	
        recoil+=aJet;
	
	nEvent.jetpt[nEvent.jetNum]  = aJet.Pt();
	nEvent.jeteta[nEvent.jetNum] = aJet.Eta();
	nEvent.jetphi[nEvent.jetNum] = aJet.Phi();
	nEvent.jetpx[nEvent.jetNum]  = aJet.Px();
	nEvent.jetpy[nEvent.jetNum]  = aJet.Py();
	nEvent.jetpz[nEvent.jetNum]  = aJet.Pz();
	nEvent.jetenergy[nEvent.jetNum] = aJet.E();
	nEvent.jetscale[nEvent.jetNum]  = jesC;
	if(isJetID) nEvent.jetID[nEvent.jetNum] = 1;
	
	nEvent.jetNum = nEvent.jetNum + 1 ;
	
        if ( jpt>30 ) {
          counters[JE].fill("... pt>30");
          nEvent.goodJetNum++;
        }
		
      }
    
    
    TLorentzVector sumOfPFJets(0,0,0,0);
    nEvent.pfJetNum=0;
    nEvent.pfJetGoodNum=0;
    for(int i =0 ; i<fTR->PFNJets;i++) // jet loop
      {
        counters[PJ].fill("All PF jets");
	if(i==jMax){cout<<"max Num was reached"<<endl; return;}
	
	float jpt = fTR->PFJPt[i];
	float jeta = fTR->PFJEta[i];
	float jphi = fTR->PFJPhi[i];
	float jpx = fTR->PFJPx[i];
	float jpy = fTR->PFJPy[i];
	float jpz = fTR->PFJPz[i];
	float jenergy = fTR->PFJE[i];
	float jesC = fTR->PFJScale[i];
	bool isJetID = IsGoodBasicPFJet(i,false);
	
	if (!isJetID) continue;
        counters[PJ].fill("... pass Jet ID");

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
	
	
	nEvent.pfJetPt[nEvent.pfJetNum] = jpt;
	nEvent.pfJetEta[nEvent.pfJetNum] = jeta;
	nEvent.pfJetPhi[nEvent.pfJetNum] = jphi;
	nEvent.pfJetNum = nEvent.pfJetNum +1;
	nEvent.pfHT += jpt;
	
	//This was added by Pablo
	if(isJetID && abs(jeta)<2.6 && jpt>15)
	  {
            counters[PJ].fill("... pass loose jet selection");
	    nEvent.pfGoodHT += jpt;
	    sumOfPFJets += aJet;
	    if ( jpt>30 ) {
              counters[PJ].fill("... pass tight jet selection");
	      nEvent.pfTightHT += jpt;
	      nEvent.pfJetGoodPt[nEvent.pfJetGoodNum]  = jpt;
	      nEvent.pfJetGoodEta[nEvent.pfJetGoodNum] = jeta;
	      nEvent.pfJetGoodPhi[nEvent.pfJetGoodNum] = jphi;
	      nEvent.pfJetGoodNum++;
	    }
	  }
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
      if ( fTR->PFMET>50 && nEvent.pfTightHT>100 ) {
        if ( fabs(nEvent.mll-91)>15 || nEvent.id1*nEvent.id2 == 2 )
          counters[EV].fill("... "+type+" + 2 jets + PFMET>50 + HT>100 + Z veto");
        else
          counters[EV].fill("... "+type+" + 2 jets + PFMET>50 + HT>100 + require Z");
      }
    }
    ////////////////////////////////////////////////////

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
    nEvent.met[1]=fTR->MuCorrMET;
    nEvent.met[2]=fTR->TCMET;
    nEvent.met[3]=fTR->MuJESCorrMET;
    nEvent.met[4]=fTR->PFMET;
    nEvent.met[5]=fTR->SumEt;
    
    nEvent.eventNum = fTR->Event;
    nEvent.runNum = fTR->Run;
    nEvent.lumi = fTR->LumiSection;
    
    nEvent.totEvents = fTR->GetEntries();
    
    float caloMETpx = fTR->RawMETpx;
    float caloMETpy = fTR->RawMETpy;
    
    float pfMETpx = fTR->PFMETpx;
    float pfMETpy = fTR->PFMETpy;
    
    TLorentzVector caloMETvector(caloMETpx,caloMETpy,0,0);
    TLorentzVector pfMETvector(pfMETpx,pfMETpy,0,0);
    
    TLorentzVector caloVector(0,0,0,0); // for constructing SumJPt from raw calomet
    TLorentzVector pfJetVector(0,0,0,0); // for constructing SumJPt from pf jets, as Pablo
    TLorentzVector pfNoCutsJetVector(0,0,0,0); // for constructing SumJPt from pfmet (unclustered), as Kostas
    
    
    
    
    if(sortedGoodLeptons[PosLepton1].type == 0 && sortedGoodLeptons[PosLepton2].type == 0)
      {
	caloVector = -caloMETvector - s1 -s2; // subtract the electrons
	pfNoCutsJetVector = -pfMETvector - s1 - s2; // remove the electrons
	
      }
    
    if(sortedGoodLeptons[PosLepton1].type == 1 && sortedGoodLeptons[PosLepton2].type == 1)
      {
	//caloVector = -caloMETvector + s1 + s2; // add the muons
	caloVector = -caloMETvector; // add the muons -> NO!
	pfNoCutsJetVector = -pfMETvector - s1 - s2; // remove the muons
      }
    
    
    // #--- different versions of JZB
    
    nEvent.dphi_sumJetVSZ[0]=caloVector.DeltaPhi(s1+s2); // DPhi between Z and SumJpt
    nEvent.sumJetPt[0]=caloVector.Pt();
    nEvent.jzb[0] = caloVector.Pt() - (s1+s2).Pt(); // calib issue of rawcalomet wrt lepton energy scale, under develop
    
    nEvent.dphi_sumJetVSZ[1] = pfNoCutsJetVector.DeltaPhi(s1+s2); 
    nEvent.sumJetPt[1] = pfNoCutsJetVector.Pt(); 
    nEvent.jzb[1] = pfNoCutsJetVector.Pt() - (s1+s2).Pt(); // to be used with pfMET
    
    
    nEvent.dphi_sumJetVSZ[2] = recoil.DeltaPhi(s1+s2);  // recoil is not yet a recoil but the sumJPt, since the leptons will be added only later (ugly)
    nEvent.sumJetPt[2] = recoil.Pt(); 
    nEvent.jzb[2] = recoil.Pt() - (s1+s2).Pt(); // to be used recoil met (recoilpt[0])
    
    
    nEvent.jzb[3] = sumOfPFJets.Pt() - (s1+s2).Pt(); // to be used recoil met (recoilpt[0])
    nEvent.sumJetPt[3] = sumOfPFJets.Pt();
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
    
    myTree->Fill();
  }
    
}

void JZBAnalysis::End(){
  fHistFile->cd();	
  /*
    fHElectronPtEta->Write();
    fHElectronIDPtEta->Write();
    fHElectronIDIsoPtEta->Write();
    for(int i =0;i<20;i++)fHMee[i]->Write();
    fHMeeDPhi->Write();
    fHMeePt->Write();
    fHMDPhiPt->Write();
    fHMZPtJ1Pt->Write();
    fHMZPtuJ1Pt->Write();
  */

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
  if ( !(fTR->MuPt[index] > 20) )       return false;
  counters[MU].fill(" ... pt > 20");
  if ( !(fabs(fTR->MuEta[index])<2.5) ) return false;
  counters[MU].fill(" ... |eta| < 2.5");

  // Quality cuts
  if ( !fTR->MuIsGMPT[index] )        return false;
  counters[MU].fill(" ... is global muon prompt tight");
  if ( !fTR->MuIsGlobalMuon[index] )  return false;
  counters[MU].fill(" ... is global muon");
  if ( !fTR->MuIsTrackerMuon[index] ) return false;
  counters[MU].fill(" ... is tracker muon");
  if ( !(fTR->MuPtE[index]/fTR->MuPt[index] < 0.1) ) return false;
  counters[MU].fill(" ... dpt/pt < 0.1");
  
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
  double hybridIso = fTR->MuRelIso03[index]*fTR->MuPt[index]/std::max(20.,fTR->MuPt[index]);
  if ( !(hybridIso < 0.15) ) return false;
  counters[MU].fill(" ... hybridIso < 0.15");

  return true;
}


const bool JZBAnalysis::IsCustomEl(const int index){

  // kinematic acceptance
  if(!(fTR->ElPt[index]>20) )return false;
  counters[EL].fill(" ... pt > 20");
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
  if(elIDWP95!=7) return false;
  counters[EL].fill(" ... passes WP95 ID");

  // Flat isolation below 20 GeV (only for synch.)
  double hybridIso = fTR->ElRelIso03[index]
    *fTR->ElPt[index]/std::max(20.,fTR->ElPt[index]);
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

  if ( !(fTR->JID_n90Hits[index] > 1) ) return false;
  counters[JE].fill(" ... n90Hits > 1");
  if ( !(fTR->JID_HPD[index] < 0.98)  ) return false;
  counters[JE].fill(" ... HPD < 0.98");

  if ( fabs(fTR->JEta[index])<2.6 ) {
    if ( !(fTR->JEMfrac[index] > 0.01)  ) return false;
  } else {
    if ( !(fTR->JEMfrac[index] > -0.9)  ) return false;
    if ( fTR->JPt[index] > 80 && !(fTR->JEMfrac[index]<1) ) return false;
  }
  counters[JE].fill(" ... pass EMfrac cut");

  return true;
}

