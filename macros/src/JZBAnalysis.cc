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
    float vjetpt[rMax];
    float vjeteta[rMax];
    float vjetphi[rMax];
    float recoilenergy[rMax];
    float recoilphi[rMax];
    float recoileta[rMax];

    float met[metMax];
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

    for(int rCounter=0;rCounter<rMax;rCounter++){
    recoilpt[rCounter]=0;
    vjetpt[rCounter]=0;
    vjeteta[rCounter]=0;
    vjetphi[rCounter]=0;
    recoilenergy[rCounter]=0;
    recoilphi[rCounter]=0;
    recoileta[rCounter]=0;
    }

    
    for(int metCounter=0;metCounter<metMax;metCounter++){
    met[metMax]=0;
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


JZBAnalysis::JZBAnalysis(TreeReader *tr, bool isData ) : UserAnalysisBase(tr), fIsData_(isData) {
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
	myTree->Branch("recoilphi",nEvent.recoilphi,"recoilphi[30]/F");
	myTree->Branch("recoileta",nEvent.recoileta,"recoileta[30]/F");
	myTree->Branch("recoilenergy",nEvent.recoilenergy,"recoilenergy[30]/F");

	myTree->Branch("vjetpt",nEvent.vjetpt,"vjetpt[30]/F");
	myTree->Branch("vjeteta",nEvent.vjeteta,"vjeteta[30]/F");
	myTree->Branch("vjetphi",nEvent.vjetphi,"vjetphi[30]/F");

	myTree->Branch("met",nEvent.met,"met[30]/F");


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
const bool JZBAnalysis::passTrigger() {
  
  
  if(  GetHLTResult("HLT_Ele15_SW_L1R") || 
       GetHLTResult("HLT_Ele15_SW_CaloEleId_L1R") ||
       GetHLTResult("HLT_Ele17_SW_CaloEleId_L1R") || 
       GetHLTResult("HLT_Ele17_SW_TightEleId_L1R") ||
       GetHLTResult("HLT_Ele17_SW_TightCaloEleId_Ele8HE_L1R_v1") ) {
    counters[EV].fill("... pass electron triggers");
    return true;
  }
  if (GetHLTResult("HLT_Mu9") || GetHLTResult("HLT_Mu15_v1")) {
    counters[EV].fill("... pass muon triggers");
    return true;
  }

  return false;
    
} 




void JZBAnalysis::Analyze(){
  
  counters[EV].fill("All events");
  if ( fIsData_ ) {
    if(!passTrigger()) return;
    counters[EV].fill("... pass all trigger requirements");
  }

  
  // #--- analysis global parameters
  double DRmax=0.4; // veto jets in a cone of DRmax close to the lepton
  nEvent.reset();
  
  // #--- trigger selection, vtx requirement
  // #--- Vertex info
  nEvent.numVtx = fTR->NVrtx;

  float rho = sqrt(fTR->PrimVtxx*fTR->PrimVtxx + fTR->PrimVtxy*fTR->PrimVtxy);
  if(fTR->PrimVtxGood)nEvent.goodVtx |=2; // save bits of vertex quality
  if(fTR->PrimVtxGood==0 && fTR->PrimVtxIsFake==0 && fTR->PrimVtxNdof>4 && fabs(fTR->PrimVtxz)<24 && rho<2)nEvent.goodVtx |=4;
  
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
    
    
    // #--- ugly part of the code starts here, needs heavy refurbishment !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    vector<TLorentzVector> Jets; // barrel+endcaps+HF
    vector<TLorentzVector> uJets;// barrel+endcap+HF uncorrected
    vector<TLorentzVector> barrelJets;// only barrel
    
    vector<TLorentzVector> pfJets; // all
    vector<TLorentzVector> pfJetGoods; // |eta|<3.0 && pfJetID
    
    // #--- construct different recoil models, initial the recoil vector will hold only the sum over the hard jets, only in the end we will add-up the lepton system
    TLorentzVector recoil(0,0,0,0); // different constructions of recoil model (under dev, need cleaning)
    TLorentzVector recoilUC(0,0,0,0); // uncorrected jets recoil
    TLorentzVector recoilID(0,0,0,0); // jets with differnt jetID (for the moment equivalent to the simple recoil 2 lines before, since the same jetID is always required)
    TLorentzVector recoilUCID(0,0,0,0); // uncorrected jets with different JetID and 
    
    TLorentzVector barrelRecoil(0,0,0,0); // same as before but only for the |eta|<1.2
    TLorentzVector barrelRecoilUC(0,0,0,0);
    TLorentzVector barrelRecoilID(0,0,0,0);
    TLorentzVector barrelRecoilUCID(0,0,0,0);
    
    TLorentzVector sumOfPFJets(0,0,0,0);

    nEvent.centraljetNum=0; // barrel jet counting
    nEvent.jetNum=0; // total jet counting
    for(int i =0 ; i<fTR->NJets;i++) // jet loop
      {
	if(i==jMax){cout<<"max Num was reached"<<endl; return;}
	
	float jpt = fTR->JPt[i];
	float jeta = fTR->JEta[i];
	float jpx = fTR->JPx[i];
	float jpy = fTR->JPy[i];
	float jpz = fTR->JPz[i];
	float jenergy = fTR->JE[i];
	float jesC = fTR->JEcorr[i];
	bool isJetID = IsCustomJet(i);
	
	if(!isJetID)continue; //consider only Jets with JetID
	
	TLorentzVector aJet(jpx,jpy,jpz,jenergy);
	
	if(aJet.DeltaR(sortedGoodLeptons[PosLepton1].p)<DRmax)continue;  // do not consider jets 2 close to the electrons
	if(aJet.DeltaR(sortedGoodLeptons[PosLepton2].p)<DRmax)continue;
	
	
	recoil+=aJet;
	recoilUC+=aJet*(1/jesC);
	
	nEvent.jetpt[nEvent.jetNum]=aJet.Pt();
	nEvent.jeteta[nEvent.jetNum]=aJet.Eta();
	nEvent.jetphi[nEvent.jetNum]=aJet.Phi();
	nEvent.jetpx[nEvent.jetNum]=aJet.Px();
	nEvent.jetpy[nEvent.jetNum]=aJet.Py();
	nEvent.jetpz[nEvent.jetNum]=aJet.Pz();
	nEvent.jetenergy[nEvent.jetNum]=aJet.E();
	nEvent.jetscale[nEvent.jetNum]=jesC;
	if(isJetID)nEvent.jetID[nEvent.jetNum]=1;
	
	nEvent.jetNum = nEvent.jetNum + 1 ;
	
	if(isJetID)
	  {
	    recoilID+=aJet;
	    recoilUCID+=aJet*(1/jesC);
	    
	  }
	
	if(abs(jeta)<1.2 && jpt>20)
	  {
	    
	    nEvent.centraljetpt[nEvent.centraljetNum]=aJet.Pt();
	    nEvent.centraljeteta[nEvent.centraljetNum]=aJet.Eta();
	    nEvent.centraljetphi[nEvent.centraljetNum]=aJet.Phi();
	    nEvent.centraljetpx[nEvent.centraljetNum]=aJet.Px();
	    nEvent.centraljetpy[nEvent.centraljetNum]=aJet.Py();
	    nEvent.centraljetpz[nEvent.centraljetNum]=aJet.Pz();
	    nEvent.centraljetenergy[nEvent.centraljetNum]=aJet.E();
	    nEvent.centraljetscale[nEvent.centraljetNum]=jesC;
	    if(isJetID)nEvent.centraljetID[nEvent.centraljetNum]=1;
	    
	    nEvent.centraljetNum = nEvent.centraljetNum + 1 ;
	    
	    barrelRecoil+=aJet;
	    barrelRecoilUC+=aJet*(1/jesC);
	    if(isJetID)
	      {
		barrelRecoilID+=aJet;
		barrelRecoilUCID+=aJet*(1/jesC);
		
	      }
	    
	  }
	
      }
    
    
    nEvent.pfJetNum=0;
    nEvent.pfJetGoodNum=0;
    for(int i =0 ; i<fTR->PFNJets;i++) // jet loop
      {
        counters[JE].fill("All PF jets");
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
	
	if (isJetID) counters[JE].fill("... pass Jet ID");

	TLorentzVector aJet(jpx,jpy,jpz,jenergy);
	
	if(aJet.DeltaR(sortedGoodLeptons[PosLepton1].p)<DRmax)continue;  // do not consider jets 2 close to the electrons
        counters[JE].fill("... pass lepton 1 veto");
	if(aJet.DeltaR(sortedGoodLeptons[PosLepton2].p)<DRmax)continue;
	counters[JE].fill("... pass lepton 2 veto");

	nEvent.pfJetPt[nEvent.pfJetNum] = jpt;
	nEvent.pfJetEta[nEvent.pfJetNum] = jeta;
	nEvent.pfJetPhi[nEvent.pfJetNum] = jphi;
	nEvent.pfJetNum = nEvent.pfJetNum +1;
	nEvent.pfHT += jpt;
	
	//This was added by Pablo
	if(isJetID && abs(jeta)<3.0 && jpt>15)
	  {
            counters[JE].fill("... pass loose jet selection");
	    nEvent.pfGoodHT += jpt;
	    sumOfPFJets += aJet;
	    if ( jpt>30 ) {
              counters[JE].fill("... pass tight jet selection");
	      nEvent.pfTightHT += jpt;
	      nEvent.pfJetGoodPt[nEvent.pfJetGoodNum]  = jpt;
	      nEvent.pfJetGoodEta[nEvent.pfJetGoodNum] = jeta;
	      nEvent.pfJetGoodPhi[nEvent.pfJetGoodNum] = jphi;
	      nEvent.pfJetGoodNum++;
	    }
	  }
	
	//This is what Kostas had
	//if(isJetID && abs(jeta)<3.0 && jpt>30)
	//  {
	//    nEvent.pfJetGoodPt[nEvent.pfJetGoodNum] = jpt;
	//    nEvent.pfJetGoodEta[nEvent.pfJetGoodNum] = jeta;
	//    nEvent.pfJetGoodPhi[nEvent.pfJetGoodNum] = jphi;
	//    nEvent.pfJetGoodNum = nEvent.pfJetGoodNum +1;
	//    nEvent.pfGoodHT += jpt;
	//  }
	
	
	
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
    
    if(recoilUC.Pt()!=0)
      {
	index=1;
	nEvent.vjetpt[index]=recoilUC.Pt();
	nEvent.vjeteta[index]=recoilUC.Eta();
	nEvent.vjetphi[index]=recoilUC.Phi();
      }
    
    if(recoilID.Pt()!=0)
      {
	index=2;
	nEvent.vjetpt[index]=recoilID.Pt();
	nEvent.vjeteta[index]=recoilID.Eta();
	nEvent.vjetphi[index]=recoilID.Phi();
      }
    
    if(recoilUCID.Pt()!=0)
      {
	index=3;
	nEvent.vjetpt[index]=recoilUCID.Pt();
	nEvent.vjeteta[index]=recoilUCID.Eta();
	nEvent.vjetphi[index]=recoilUCID.Phi();
      }
    
    if(barrelRecoil.Pt()!=0)
      {
	index=4;
	nEvent.vjetpt[index]=barrelRecoil.Pt();
	nEvent.vjeteta[index]=barrelRecoil.Eta();
	nEvent.vjetphi[index]=barrelRecoil.Phi();
      }
    
    if(barrelRecoilUC.Pt()!=0)
      {
	index=5;
	nEvent.vjetpt[index]=barrelRecoilUC.Pt();
	nEvent.vjeteta[index]=barrelRecoilUC.Eta();
	nEvent.vjetphi[index]=barrelRecoilUC.Phi();
      }
    
    if(barrelRecoilID.Pt()!=0)
      {
	index=6;
	nEvent.vjetpt[index]=barrelRecoilID.Pt();
	nEvent.vjeteta[index]=barrelRecoilID.Eta();
	nEvent.vjetphi[index]=barrelRecoilID.Phi();
      }
    
    if(barrelRecoilUCID.Pt()!=0)
      {
	index=7;
	nEvent.vjetpt[index]=barrelRecoilUCID.Pt();
	nEvent.vjeteta[index]=barrelRecoilUCID.Eta();
	nEvent.vjetphi[index]=barrelRecoilUCID.Phi();
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
	caloVector = -caloMETvector + s1 + s2; // add the muons
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
    
    
    // ----------------------------------------
    recoil+=s1+s2;   // now add also the leptons to the recoil! to form the complete recoil model
    recoilUC+=s1+s2;
    recoilID+=s1+s2;
    recoilUCID+=s1+s2;
    
    barrelRecoil+=s1+s2;
    barrelRecoilUC+=s1+s2;
    barrelRecoilID+=s1+s2;
    barrelRecoilUCID+=s1+s2;
    
    if(recoil.Pt()!=0)
      {
	index=0;
	nEvent.recoilpt[index]=recoil.Pt();
	nEvent.recoileta[index]=recoil.Eta();
	nEvent.recoilphi[index]=recoil.Phi();
      }
    
    if(recoilUC.Pt()!=0)
      {
	index=1;
	nEvent.recoilpt[index]=recoilUC.Pt();
	nEvent.recoileta[index]=recoilUC.Eta();
	nEvent.recoilphi[index]=recoilUC.Phi();
      }
    
    if(recoilID.Pt()!=0)
      {
	index=2;
	nEvent.recoilpt[index]=recoilID.Pt();
	nEvent.recoileta[index]=recoilID.Eta();
	nEvent.recoilphi[index]=recoilID.Phi();
      }
    
    if(recoilUCID.Pt()!=0)
      {
	index=3;
	nEvent.recoilpt[index]=recoilUCID.Pt();
	nEvent.recoileta[index]=recoilUCID.Eta();
	nEvent.recoilphi[index]=recoilUCID.Phi();
      }
    
    if(barrelRecoil.Pt()!=0)
      {
	index=4;
	nEvent.recoilpt[index]=barrelRecoil.Pt();
	nEvent.recoileta[index]=barrelRecoil.Eta();
	nEvent.recoilphi[index]=barrelRecoil.Phi();
      }
    
    if(barrelRecoilUC.Pt()!=0)
      {
	index=5;
	nEvent.recoilpt[index]=barrelRecoilUC.Pt();
	nEvent.recoileta[index]=barrelRecoilUC.Eta();
	nEvent.recoilphi[index]=barrelRecoilUC.Phi();
      }
    
    if(barrelRecoilID.Pt()!=0)
      {
	index=6;
	nEvent.recoilpt[index]=barrelRecoilID.Pt();
	nEvent.recoileta[index]=barrelRecoilID.Eta();
	nEvent.recoilphi[index]=barrelRecoilID.Phi();
      }
    
    if(barrelRecoilUCID.Pt()!=0)
      {
	index=7;
	nEvent.recoilpt[index]=barrelRecoilUCID.Pt();
	nEvent.recoileta[index]=barrelRecoilUCID.Eta();
	nEvent.recoilphi[index]=barrelRecoilUCID.Phi();
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
  if ( !fTR->MuIsGMPT[index] )        return false;
  counters[MU].fill("... is global muon prompt tight");
  if ( !fTR->MuIsGlobalMuon[index] )  return false;
  counters[MU].fill("... is global muon");
  if ( !fTR->MuIsTrackerMuon[index] ) return false;
  counters[MU].fill("... is tracker muon");
  
  if ( !(fTR->MuPt[index] > 10) )       return false;
  counters[MU].fill("... pt > 10");
  if ( !(fabs(fTR->MuEta[index])<2.5) ) return false;
  counters[MU].fill("... |eta|<2.5");

  if ( !(fTR->MuNTkHits[index] >= 11) )     return false;
  counters[MU].fill("... nTkHits >= 11");
  if ( !(fTR->MuNPxHits[index] > 0) )       return false;
  counters[MU].fill("... nPxHits > 0");
  if ( !(fTR->MuNMatches[index] > 1) )      return false;
  counters[MU].fill("... nMatches > 1");

  if ( !(fabs(fTR->MuD0PV[index]) < 0.02) ) return false;
  counters[MU].fill("... D0(pv) < 0.02");
  if ( !(fabs(fTR->MuDzPV[index]) < 1.0 ) ) return false;
  counters[MU].fill("... DZ(pv) < 1.0");

  // Flat isolation below 20 GeV (only for synch.)
  double hybridIso = fTR->MuRelIso03[index]*fTR->MuPt[index]/std::max(20.,fTR->MuPt[index]);
  if ( !(hybridIso < 0.15) ) return false;
  counters[MU].fill("... hybridIso < 0.15");

  return true;
}


const bool JZBAnalysis::IsCustomEl(const int index){

  // kinematic acceptance
  if(!(fTR->ElPt[index]>10) )return false;
  counters[EL].fill("... pt>10");
  if(!(fabs(fTR->ElEta[index]) < 2.5) ) return false;
  counters[EL].fill("... |eta| < 2.5");
  if ( !(fTR->ElNumberOfMissingInnerHits[index] <= 1 ) ) return false;
  counters[EL].fill("... missing inner hits <= 1");
  if ( !(fabs(fTR->ElD0PV[index]) < 0.04) ) return false;
  counters[EL].fill("... D0(pv) < 0.04");
  if ( !(fabs(fTR->ElDzPV[index]) < 1.0 ) ) return false;
  counters[EL].fill("... DZ(pv) < 1.0");

  // Electron ID
  int elIDWP95 = fTR->ElIDsimpleWP95relIso[index];
  if(elIDWP95!=7) return false;
  counters[EL].fill("... passes WP95 ID");

  // Flat isolation below 20 GeV (only for synch.)
  double hybridIso = fTR->ElRelIso03[index]
    *fTR->ElPt[index]/std::max(20.,fTR->ElPt[index]);
  if ( !(hybridIso < 0.15) ) return false;  
  counters[EL].fill("... hybridIso < 0.15");

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
        // Basic Jet cleaning and ID cuts
//        if(fTR->JPt[index] < 30) return false;
//        if(fabs(fTR->JEta[index]) > 3.0) return false;

        if(fTR->JEt[index] - fTR->JPt[index] < -0.0001 ) return false;
        if(fTR->JID_n90Hits[index] < 2 ) return false;
        if(fTR->JID_HPD[index] > 0.98 ) return false;
        if(fTR->JID_RBX[index] > 0.95 ) return false;
        if(fTR->JEMfrac[index] > 1. ) return false;
        if(fTR->JEMfrac[index] < 0.01 ) return false;
        // Have a linearly decreasing cut value in the transition region where
        // the jet cones intersects with the tracker acceptance, i.e. between
        // eta 1.9 and 2.9
        const double chmin = 0.05;
        double temp = chmin;
        if(fabs(fTR->JEta[index]) > 1.9) temp = chmin * (1. - fabs(fTR->JEta[index]) + 1.9);
        if(fabs(fTR->JEta[index]) > 2.9) temp = 0.;
        if( fTR->JChfrac[index] < temp && fabs(fTR->JEta[index]) < 2.9) return false;
        return true;
}

