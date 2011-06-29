/*----------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------
------------------- Efficiency Calculator                ---------------------------------
-------------------                                      ---------------------------------
------------------- Author: pablom@cern.ch               ---------------------------------
------------------- Date: June 2011                      ---------------------------------
------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------*/

#include "helper/Utilities.hh"
#include "RunEfficiency.hh"
#include "TF1.h"
#include <time.h>
#include <TRandom.h>
#include "TF1.h"


using namespace std;


#define jMax 30  // do not touch this
#define metMax 30
#define rMax 30
#define nAnalysis 3


//_______________________t_nanoEvent class____________________________________________________
//____________________________________________________________________________________________
//Class t_nanoEvent is inspired in the JZB analysis but with a completely different definition
//____________________________________________________________________________________________

class t_nanoEvent {
public:
  
  t_nanoEvent();
  void reset();
  
  //Di-lepton system for both pf and reco 
  float mll; 
  float pfmll;
  float pt;
  float pfpt;
  float phi;
  float pfphi;
  float pfeta;
  
  //Leading leptons
  float pt1; 
  float pt2;
  float ptn;
  float pfpt1; 
  float pfpt2;
  float pfptn;
  float phi1;
  float phi2;
  float phin;
  float pfphi1;
  float pfphi2;
  float pfphin;
  float eta1;
  float eta2;
  float etan;
  float pfeta1;
  float pfeta2;
  float pfetan;
  float d0bs1;
  float dzbs1;
  float d0pv1;
  float dzpv1;
  float d0bs2;
  float dzbs2;
  float d0pv2;
  float dzpv2;
  float d0bsn;
  float dzbsn;
  float d0pvn;
  float dzpvn;
  float dr1[3];
  float dr2[3];
  float drn[3];
  int ch1;
  int ch2;
  int chn;
  int pfch1;
  int pfch2;
  int pfchn;
  int id1;
  int id2;
  int idn;
  
  int tagReco1[3];
  int probeReco1[3];
  int pprobeReco1[3];
  int tagReco2[3];
  int probeReco2[3];
  int pprobeReco2[3];
  
  int tagIso1[3];
  int probeIso1[3];
  int pprobeIso1[3];
  int tagIso2[3];
  int probeIso2[3];
  int pprobeIso2[3];
  
  int tagID1[3];
  int probeID1[3];
  int pprobeID1[3];
  int tagID2[3];
  int probeID2[3];
  int pprobeID2[3];

  int tag1[3];
  int probe1[3];
  int pprobe1[3];
  int tag2[3];
  int probe2[3];
  int pprobe2[3];

  float drl;
  float pfdrl;


  //MC information
  float genPt1;
  float genPt2;
  float genPtN;
  int genId1;
  int genId2;
  int genMID1;
  int genMID2;
  float genEta1;
  float genEta2;
  float genEtaN;
  float genPhi1;
  float genPhi2;
  int genCh1;
  int genCh2;
  float genMET;
  float genZPt;    
  float genMll;
  int   genNjets;
  int   genNleptons;
  float genDRN;
 
  //Jets 
  int pfJetGoodNum[3];
  int pfJetGoodNumID[3];
  float pfMET[3];
  float pfHT[3];
  float pfGoodHT[3];
  float pfTightHT[3];
  
  int eventNum;
  int runNum;
  int lumi;
  int goodVtx;
  int numVtx;
  
  int passedee_triggers[3];
  int passedmm_triggers[3];
  int passedem_triggers[3];
 
  
};

//Constructor
t_nanoEvent::t_nanoEvent(){};

void t_nanoEvent::reset() {
  
  mll = 0; pfmll = 0; pt = 0; pfpt = 0; phi = 0; pfphi = 0; pfeta = 0;
  pt1 = 0; phi1 = 0; eta1 = 0; pfpt1 = 0; pfphi1 = 0; pfeta1 = 0;
  pt2 = 0; phi2 = 0; eta2 = 0; pfpt2 = 0; pfphi2 = 0; pfeta2 = 0;
  ptn = 0; phin = 0; etan = 0; pfptn = 0; pfphin = 0; pfetan = 0;
  d0bs1 = 0; dzbs1 = 0; d0pv1 = 0; dzpv1 = 0; d0bs2 = 0; dzbs2 = 0; d0pv2 = 0; dzpv2 = 0; 
  d0bsn = 0; dzbsn = 0; d0pvn = 0; dzpvn = 0;
  ch1 = 0; ch2 = 0; chn = 0; id1 = -1; id2 = -1; idn = -1;
  for(int n = 0; n < nAnalysis; n++) {
    tag1[n] = 0; probe1[n] = 0; pprobe1[n] = 0;
    tag2[n] = 0; probe2[n] = 0; pprobe2[n] = 0;
  
    tagReco1[n] = 0; probeReco1[n] = 0; pprobeReco1[n] = 0;
    tagReco2[n] = 0; probeReco2[n] = 0; pprobeReco2[n] = 0;
  
    tagIso1[n] = 0; probeIso1[n] = 0; pprobeIso1[n] = 0;
    tagIso2[n] = 0; probeIso2[n] = 0; pprobeIso2[n] = 0;
    
    tagID1[n] = 0; probeID1[n] = 0; pprobeID1[n] = 0;
    tagID2[n] = 0; probeID2[n] = 0; pprobeID2[n] = 0;

    dr1[n] = 0; dr2[n] = 0; drn[n] = 0; 
  
    pfJetGoodNum[n] = 0; pfJetGoodNumID[n] = 0;
    pfHT[n] = 0; pfGoodHT[n] = 0; pfTightHT[n] = 0;

    passedee_triggers[n] = 0; passedmm_triggers[n] = 0; passedem_triggers[n] = 0;

    pfMET[n]=0;
  }

  drl = 0; pfdrl = 0;

  genPt1 = 0; genId1 = 0; genMID1 = 0; genEta1 = 0; genPhi1 = 0; genCh1 = 0;
  genPt2 = 0; genId2 = 0; genMID2 = 0; genEta2 = 0; genPhi2 = 0; genCh2 = 0;
  genPtN = 0; genEtaN = 0;
  genMET = 0; genZPt = 0; genMll = 0; genNjets = 0; genNleptons = 0;
  genDRN = 0;

  eventNum=0; runNum=0; lumi=0; goodVtx=0; numVtx=0;


  

}
//_____________________End of t_nanoEvent class_______________________________________
//____________________________________________________________________________________



TTree *t_myTree;
TTree *t_myInfo;

t_nanoEvent t_nEvent;


//____________________Function used by sort____________________________________________
bool momentumComparator(t_lepton i, t_lepton j) { return (i.p.Pt()>j.p.Pt()); }


//____________________Start of RunEfficiency___________________________________________
//_____________________________________________________________________________________
//_____________________________________________________________________________________
RunEfficiency::RunEfficiency(TreeReader *tr, std::string dataType, bool fullCleaning) : 
  UserAnalysisBase(tr), fDataType_(dataType), fFullCleaning_(fullCleaning) {
  //	Util::SetStyle();	
  //	setTDRStyle();	
}

//_____________________________________________________________________________________
RunEfficiency::~RunEfficiency(){
}


//_____________________________________________________________________________________
void RunEfficiency::setInfo() {
 
  TH1::AddDirectory(kFALSE);
  // Taken from JZB analysis
  t_myInfo = new TTree("info","info/S");
  TString *user = new TString();
  TString *timestamp = new TString();
  TString *cmsdir = new TString();
  t_myInfo->Branch("user",&user,16000,0);
  t_myInfo->Branch("timestamp",&timestamp,16000,0);
  t_myInfo->Branch("cmsdir",&cmsdir,16000,0);
  
  char usertext[255];
  FILE *usernamefile;
  usernamefile = popen("whoami", "r");
  fgets(usertext, sizeof(usertext), usernamefile);
  pclose(usernamefile);
  char scmsdir[1000];
  getcwd(scmsdir,1000);
  *cmsdir=scmsdir;
  *user=usertext;
  time_t rawtime;
  struct tm * timeinfo;
  time (&rawtime );
  *timestamp=ctime(&rawtime);
  
  t_myInfo->Fill();
  t_myInfo->Write();
  
}


//__________________________________________________________________________
void RunEfficiency::Begin(){

  // Define the output file of histograms
  fHistFile = new TFile(outputFileName_.c_str(), "RECREATE");
  
  setInfo();
  
  t_myTree = new TTree("events","events");
  t_myTree->Branch("mll",&t_nEvent.mll,"mll/F");
  t_myTree->Branch("pt",&t_nEvent.pt,"pt/F");
  t_myTree->Branch("phi",&t_nEvent.phi,"phi/F");
  t_myTree->Branch("pt1",&t_nEvent.pt1,"pt1/F");
  t_myTree->Branch("pt2",&t_nEvent.pt2,"pt2/F");
  t_myTree->Branch("ptn",&t_nEvent.ptn,"ptn/F");
  t_myTree->Branch("dr1",t_nEvent.dr1,"dr1[3]/F");
  t_myTree->Branch("dr2",t_nEvent.dr2,"dr2[3]/F");
  t_myTree->Branch("drn",t_nEvent.drn,"drn[3]/F");
  t_myTree->Branch("d0bs1",&t_nEvent.d0bs1,"d0bs1/F");
  t_myTree->Branch("dzbs1",&t_nEvent.dzbs1,"dzbs1/F");
  t_myTree->Branch("d0pv1",&t_nEvent.d0pv1,"d0pv1/F");
  t_myTree->Branch("dzpv1",&t_nEvent.dzpv1,"dzpv1/F");
  t_myTree->Branch("d0bs2",&t_nEvent.d0bs2,"d0bs2/F");
  t_myTree->Branch("dzbs2",&t_nEvent.dzbs2,"dzbs2/F");
  t_myTree->Branch("d0pv2",&t_nEvent.d0pv2,"d0pv2/F");
  t_myTree->Branch("dzpv2",&t_nEvent.dzpv2,"dzpv2/F");
  t_myTree->Branch("d0bsn",&t_nEvent.d0bsn,"d0bsn/F");
  t_myTree->Branch("dzbsn",&t_nEvent.dzbsn,"dzbsn/F");
  t_myTree->Branch("d0pvn",&t_nEvent.d0pvn,"d0pvn/F");
  t_myTree->Branch("dzpvn",&t_nEvent.dzpvn,"dzpvn/F");
  t_myTree->Branch("genPt1",&t_nEvent.genPt1,"genPt1/F");
  t_myTree->Branch("genPt2",&t_nEvent.genPt2,"genPt2/F");
  t_myTree->Branch("genPtN",&t_nEvent.genPtN,"genPtN/F");
  t_myTree->Branch("genEta1",&t_nEvent.genEta1,"genEta1/F");
  t_myTree->Branch("genEta2",&t_nEvent.genEta2,"genEta2/F");
  t_myTree->Branch("genEtaN",&t_nEvent.genEtaN,"genEtaN/F");
  t_myTree->Branch("genPhi1",&t_nEvent.genPhi1,"genPhi1/F");
  t_myTree->Branch("genPhi2",&t_nEvent.genPhi2,"genPhi2/F"); 
  t_myTree->Branch("genId1",&t_nEvent.genId1,"genId1/I");
  t_myTree->Branch("genId2",&t_nEvent.genId2,"genId2/I");
  t_myTree->Branch("genMID1",&t_nEvent.genMID1,"genMID1/I");
  t_myTree->Branch("genMID2",&t_nEvent.genMID2,"genMID2/I");
  t_myTree->Branch("genChD1",&t_nEvent.genCh1,"genCh1/I");
  t_myTree->Branch("genChD2",&t_nEvent.genCh2,"genCh2/I");
  t_myTree->Branch("genMET",&t_nEvent.genMET,"genMET/F");
  t_myTree->Branch("genZPt",&t_nEvent.genZPt,"genZPt/F");
  t_myTree->Branch("genMll",&t_nEvent.genMll,"genMll/F");
  t_myTree->Branch("genNjets",&t_nEvent.genNjets,"genNjets/I");
  t_myTree->Branch("genNleptons",&t_nEvent.genNleptons,"genNleptons/I");
  t_myTree->Branch("genDRN",&t_nEvent.genDRN,"genDRN/F");
  t_myTree->Branch("eta1",&t_nEvent.eta1,"eta1/F");
  t_myTree->Branch("eta2",&t_nEvent.eta2,"eta2/F");
  t_myTree->Branch("etan",&t_nEvent.etan,"etan/F");
  t_myTree->Branch("phi1",&t_nEvent.phi1,"phi1/F");
  t_myTree->Branch("phi2",&t_nEvent.phi2,"phi2/F");
  t_myTree->Branch("phin",&t_nEvent.phin,"phin/F");
  t_myTree->Branch("id1",&t_nEvent.id1,"id1/I");
  t_myTree->Branch("id2",&t_nEvent.id2,"id2/I");
  t_myTree->Branch("idn",&t_nEvent.idn,"idn/I");
  t_myTree->Branch("ch1",&t_nEvent.ch1,"ch1/I");
  t_myTree->Branch("ch2",&t_nEvent.ch2,"ch2/I");
  t_myTree->Branch("chn",&t_nEvent.chn,"chn/I");
  t_myTree->Branch("drl", &t_nEvent.drl,"drl/F");
  t_myTree->Branch("pfdrl", &t_nEvent.pfdrl,"pfdrl/F");
  t_myTree->Branch("pfMET",&t_nEvent.pfMET,"pfMET[3]/F");
  t_myTree->Branch("pfHT",t_nEvent.pfHT,"pfHT[3]/F");
  t_myTree->Branch("pfGoodHT", t_nEvent.pfGoodHT,"pfGoodHT[3]/F");
  t_myTree->Branch("pfTightHT", t_nEvent.pfTightHT,"pfTightHT[3]/F");
  t_myTree->Branch("eventNum",&t_nEvent.eventNum,"eventNum/I");
  t_myTree->Branch("runNum",&t_nEvent.runNum,"runNum/I");
  t_myTree->Branch("lumi",&t_nEvent.lumi,"lumi/I");
  t_myTree->Branch("goodVtx",&t_nEvent.goodVtx,"goodVtx/I");
  t_myTree->Branch("numVtx",&t_nEvent.numVtx,"numVtx/I");
  t_myTree->Branch("pfJetGoodNum",t_nEvent.pfJetGoodNum,"pfJetGoodNum[3]/I");
  t_myTree->Branch("pfJetGoodNumID",t_nEvent.pfJetGoodNumID,"pfJetGoodNumID[3]/I");
  t_myTree->Branch("pfpt",&t_nEvent.pfpt,"pfpt/F");
  t_myTree->Branch("pfphi",&t_nEvent.pfphi,"pfphi/F");
  t_myTree->Branch("pfeta",&t_nEvent.pfeta,"pfeta/F");
  t_myTree->Branch("pfpt1",&t_nEvent.pfpt1,"pfpt1/F");
  t_myTree->Branch("pfpt2",&t_nEvent.pfpt2,"pfpt2/F");
  t_myTree->Branch("pfptn",&t_nEvent.pfptn,"pfptn/F");
  t_myTree->Branch("pfphi1",&t_nEvent.pfphi1,"pfphi1/F");
  t_myTree->Branch("pfphi2",&t_nEvent.pfphi2,"pfphi2/F");
  t_myTree->Branch("pfphin",&t_nEvent.pfphin,"pfphin/F");
  t_myTree->Branch("pfeta1",&t_nEvent.pfeta1,"pfeta1/F");
  t_myTree->Branch("pfeta2",&t_nEvent.pfeta2,"pfeta2/F");
  t_myTree->Branch("pfetan",&t_nEvent.pfetan,"pfetan/F");
  t_myTree->Branch("pfmll",&t_nEvent.pfmll,"pfmll/F");
  t_myTree->Branch("passedee_triggers", &t_nEvent.passedee_triggers,"passedee_triggers[3]/I");
  t_myTree->Branch("passedmm_triggers", &t_nEvent.passedmm_triggers,"passedmm_triggers[3]/I");
  t_myTree->Branch("passedem_triggers", &t_nEvent.passedem_triggers,"passedem_triggers[3]/I");
  t_myTree->Branch("tag1", t_nEvent.tag1,"tag1[3]/I");
  t_myTree->Branch("probe1", t_nEvent.probe1,"probe1[3]/I");
  t_myTree->Branch("pprobe1", t_nEvent.pprobe1,"pprobe1[3]/I");
  t_myTree->Branch("tag2", t_nEvent.tag2,"tag2[3]/I");
  t_myTree->Branch("probe2", t_nEvent.probe2,"probe2[3]/I");
  t_myTree->Branch("pprobe2", t_nEvent.pprobe2,"pprobe2[3]/I");
  t_myTree->Branch("tagReco1", t_nEvent.tagReco1,"tagReco1[3]/I");
  t_myTree->Branch("probeReco1", t_nEvent.probeReco1,"probeReco1[3]/I");
  t_myTree->Branch("pprobeReco1", t_nEvent.pprobeReco1,"pprobeReco1[3]/I");
  t_myTree->Branch("tagReco2", t_nEvent.tagReco2,"tagReco2[3]/I");
  t_myTree->Branch("probeReco2", t_nEvent.probeReco2,"probeReco2[3]/I");
  t_myTree->Branch("pprobeReco2", t_nEvent.pprobeReco2,"pprobeReco2[3]/I");
  t_myTree->Branch("tagIso1", t_nEvent.tagIso1,"tagIso1[3]/I");
  t_myTree->Branch("probeIso1", t_nEvent.probeIso1,"probeIso1[3]/I");
  t_myTree->Branch("pprobeIso1", t_nEvent.pprobeIso1,"pprobeIso1[3]/I");
  t_myTree->Branch("tagIso2", t_nEvent.tagIso2,"tagIso2[3]/I");
  t_myTree->Branch("probeIso2", t_nEvent.probeIso2,"probeIso2[3]/I");
  t_myTree->Branch("pprobeIso2", t_nEvent.pprobeIso2,"pprobeIso2[3]/I");
  t_myTree->Branch("tagID1", t_nEvent.tagID1,"tagID1[3]/I");
  t_myTree->Branch("probeID1", t_nEvent.probeID1,"probeID1[3]/I");
  t_myTree->Branch("pprobeID1", t_nEvent.pprobeID1,"pprobeID1[3]/I");
  t_myTree->Branch("tagID2", t_nEvent.tagID2,"tagID2[3]/I");
  t_myTree->Branch("probeID2", t_nEvent.probeID2,"probeID2[3]/I");
  t_myTree->Branch("pprobeID2", t_nEvent.pprobeID2,"pprobeID2[3]/I");
  
}



//________________________________________________________________________________
vector<t_lepton> RunEfficiency::sortLeptonsByPt(vector<t_lepton>& leptons) {
  
  vector<t_lepton> theLep = leptons;
  sort (theLep.begin(), theLep.end(), momentumComparator);
  return theLep;  
  
}


//________________________________________________________________________________
const bool RunEfficiency::passElTriggers(int iAnalysis) {
  if(iAnalysis==1) {
	if ( GetHLTResult("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v1") )        return true;
	if ( GetHLTResult("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v2") )        return true;
	if ( GetHLTResult("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v3") )        return true;
	if ( GetHLTResult("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v4") )        return true;
	if ( GetHLTResult("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v5") )        return true;
   }
  return false;

}

//________________________________________________________________________________
const bool RunEfficiency::passMuTriggers(int iAnalysis) {
  if(iAnalysis==1) {
	if ( GetHLTResult("HLT_DoubleMu6_v1") )        return true;
	if ( GetHLTResult("HLT_DoubleMu6_v2") )        return true;
	if ( GetHLTResult("HLT_DoubleMu6_v3") )        return true;
	if ( GetHLTResult("HLT_DoubleMu7_v1") )        return true;
	if ( GetHLTResult("HLT_DoubleMu7_v2") )        return true;
	if ( GetHLTResult("HLT_DoubleMu7_v3") )        return true;
  }
  return false;
} 

//______________________________________________________________________________
const bool RunEfficiency::passEMuTriggers(int iAnalysis) {
  if(iAnalysis==1) {
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
  }
  return false;
}


//______________________________________________________________________________
void RunEfficiency::Analyze() {

  // #--- analysis global parameters
  double DRmax=0.4; // veto jets in a cone of DRmax close to the lepton


  t_nEvent.reset();
  
  bool isMC = (fDataType_ == "mc");
  bool isZ = false;

  //Triggers
  if ( !isMC ) {
    for(int n = 0; n < nAnalysis; n++) {
	if(passElTriggers(n)) t_nEvent.passedee_triggers[n] = 1;
	if(passMuTriggers(n)) t_nEvent.passedmm_triggers[n] = 1;
	if(passEMuTriggers(n)) t_nEvent.passedem_triggers[n] = 1;
    }
  } else {
    isZ = GeneratorInfo();
    if(!isZ ) return;
  }

  
  // #--- Vertex info
  // Fill generic information
  t_nEvent.eventNum  = fTR->Event;
  t_nEvent.runNum    = fTR->Run;
  t_nEvent.lumi      = fTR->LumiSection;
  t_nEvent.numVtx = fTR->NVrtx;
  
  float rho = sqrt(fTR->PrimVtxx*fTR->PrimVtxx + fTR->PrimVtxy*fTR->PrimVtxy);
  if(fTR->PrimVtxGood) t_nEvent.goodVtx |=2; // save bits of vertex quality
  if (   fTR->PrimVtxGood==0 && fTR->PrimVtxIsFake==0 
	 && fTR->PrimVtxNdof>4  && fabs(fTR->PrimVtxz)<24 && rho<2)
    t_nEvent.goodVtx |=4;
  
  // Good event requirement: essentially vertex requirements
  if ( !IsGoodEvent() ) {
    if( isZ ) t_myTree->Fill();
    return;
  }
  
  vector<t_lepton> leptons;
  TLorentzVector genZvector; // To store the true Z vector
  
  // #--- muon loop
  for(int muIndex=0;muIndex<fTR->NMus;muIndex++) {
    
    int indexOfAssociatedPFMuon = getPFMuIndex(muIndex);
    float px= fTR->MuPx[muIndex];
    float py= fTR->MuPy[muIndex];
    float pz= fTR->MuPz[muIndex];
    float energy = fTR->MuE[muIndex];
    
    float pfpx = 0;
    float pfpy = 0;
    float pfpz = 0;
    float pfenergy = 0 ;
    if(indexOfAssociatedPFMuon != -1) {
      pfpx= fTR->PfMu3Px[indexOfAssociatedPFMuon];
      pfpy= fTR->PfMu3Py[indexOfAssociatedPFMuon];
      pfpz= fTR->PfMu3Pz[indexOfAssociatedPFMuon];
      pfenergy =  fTR->PfMu3E[indexOfAssociatedPFMuon];
    }

    TLorentzVector tmpVector(px,py,pz,energy);
    TLorentzVector tmpPfVector(pfpx,pfpy,pfpz,pfenergy);
    
    int tmpCharge = fTR->MuCharge[muIndex];
    t_lepton tmpLepton;
    tmpLepton.p = tmpVector;
    tmpLepton.pfp = tmpPfVector;
    tmpLepton.charge = tmpCharge;
    tmpLepton.index = muIndex;
    tmpLepton.type = 1;
    tmpLepton.genPt = 0.;
    tmpLepton.d0BS = fTR->MuD0BS[muIndex];
    tmpLepton.dzBS = fTR->MuDzBS[muIndex];
    tmpLepton.d0PV = fTR->MuD0PV[muIndex];
    tmpLepton.dzPV = fTR->MuD0PV[muIndex];
    tmpLepton.pfindex = indexOfAssociatedPFMuon;
    for(int n = 0; n < nAnalysis; n++) { 
      tmpLepton.tagReco[n] = 0;
      tmpLepton.probeReco[n] = 0;
      tmpLepton.pprobeReco[n] = 0;
      tmpLepton.tagIso[n] = 0;
      tmpLepton.probeIso[n] = 0;
      tmpLepton.pprobeIso[n] = 0;
      tmpLepton.tagID[n] = 0;
      tmpLepton.probeID[n] = 0;
      tmpLepton.pprobeID[n] = 0;

      if(MuPassingRecoTag(n, muIndex)) tmpLepton.tagReco[n] = 1;
      if(MuPassingRecoProbe(n, muIndex)) tmpLepton.probeReco[n] = 1;
      if(MuPassingRecoPProbe(n, muIndex, indexOfAssociatedPFMuon)) tmpLepton.pprobeReco[n] = 1;
      
      if(MuPassingIsoTag(n, muIndex)) tmpLepton.tagIso[n] = 1;
      if(MuPassingIsoProbe(n, muIndex)) tmpLepton.probeIso[n] = 1;
      if(MuPassingIsoPProbe(n, muIndex, indexOfAssociatedPFMuon)) tmpLepton.pprobeIso[n] = 1;
      
      if(MuPassingIDTag(n, muIndex)) tmpLepton.tagID[n] = 1;
      if(MuPassingIDProbe(n, muIndex)) tmpLepton.probeID[n] = 1;
      if(MuPassingIDPProbe(n, muIndex, indexOfAssociatedPFMuon)) tmpLepton.pprobeID[n] = 1;
    }
    leptons.push_back(tmpLepton);
    
  }
  
  
  // #--- electron loop
  for(int elIndex=0;elIndex<fTR->NEles;elIndex++) {
    
    int indexOfAssociatedPFEl = getPFElIndex(elIndex);
    float px= fTR->ElPx[elIndex];
    float py= fTR->ElPy[elIndex];
    float pz= fTR->ElPz[elIndex];
    float energy =  fTR->ElE[elIndex];
    
    float pfpx = 0;
    float pfpy = 0;
    float pfpz = 0;
    float pfenergy = 0 ;
    if(indexOfAssociatedPFEl != -1) {
      pfpx= fTR->PfEl3Px[indexOfAssociatedPFEl];
      pfpy= fTR->PfEl3Py[indexOfAssociatedPFEl];
      pfpz= fTR->PfEl3Pz[indexOfAssociatedPFEl];
      pfenergy =  fTR->PfEl3E[indexOfAssociatedPFEl];
    }

    
    TLorentzVector tmpVector(px,py,pz,energy);
    TLorentzVector tmpPfVector(pfpx,pfpy,pfpz,pfenergy);
    int tmpCharge=fTR->ElCharge[elIndex];
    t_lepton tmpLepton;
    tmpLepton.pfindex = indexOfAssociatedPFEl;
    tmpLepton.p = tmpVector;
    tmpLepton.pfp = tmpPfVector;
    tmpLepton.charge = tmpCharge;
    tmpLepton.index = elIndex;
    tmpLepton.d0BS = fTR->ElD0BS[elIndex];
    tmpLepton.dzBS = fTR->ElDzBS[elIndex];
    tmpLepton.d0PV = fTR->ElD0PV[elIndex];
    tmpLepton.dzPV = fTR->ElD0PV[elIndex];
    tmpLepton.type = 0;
    tmpLepton.genPt = 0.;

    for(int n = 0; n < nAnalysis; n++) {
      tmpLepton.tagReco[n] = 0;
      tmpLepton.probeReco[n] = 0;
      tmpLepton.pprobeReco[n] = 0;
      tmpLepton.tagIso[n] = 0;
      tmpLepton.probeIso[n] = 0;
      tmpLepton.pprobeIso[n] = 0;
      tmpLepton.tagID[n] = 0;
      tmpLepton.probeID[n] = 0;
      tmpLepton.pprobeID[n] = 0;

      if(ElPassingRecoTag(n, elIndex)) tmpLepton.tagReco[n] = 1;
      if(ElPassingRecoProbe(n, elIndex)) tmpLepton.probeReco[n] = 1;
      if(ElPassingRecoPProbe(n, elIndex)) tmpLepton.pprobeReco[n] = 1;

      if(ElPassingIsoTag(n, elIndex)) tmpLepton.tagIso[n] = 1;
      if(ElPassingIsoProbe(n, elIndex)) tmpLepton.probeIso[n] = 1;
      if(ElPassingIsoPProbe(n, elIndex)) tmpLepton.pprobeIso[n] = 1;

      if(ElPassingIDTag(n, elIndex)) tmpLepton.tagID[n] = 1;
      if(ElPassingIDProbe(n, elIndex)) tmpLepton.probeID[n] = 1;
      if(ElPassingIDPProbe(n, elIndex,indexOfAssociatedPFEl)) tmpLepton.pprobeID[n] = 1;
    }
    leptons.push_back(tmpLepton);
  
  }
  

  vector<t_lepton> sortedGoodLeptons = sortLeptonsByPt(leptons);

  if(sortedGoodLeptons.size() > 1) {
    
    int PosLepton1 = 0;
    int PosLepton2 = 1;
    int negativePosition = 0;
    // Check for OS combination
    for(; PosLepton2 < sortedGoodLeptons.size(); PosLepton2++) {
      if(sortedGoodLeptons[0].charge*sortedGoodLeptons[PosLepton2].charge<0) break;
    }
    
    if(PosLepton2 == sortedGoodLeptons.size()) {
      if( isZ ) t_myTree->Fill();
      return;
    }
    
    // Acceptance
    if(sortedGoodLeptons[PosLepton1].p.Pt() > 5 && sortedGoodLeptons[PosLepton2].p.Pt() > 5 && 
       fabs(sortedGoodLeptons[PosLepton1].p.Eta())< 2.4 && fabs(sortedGoodLeptons[PosLepton2].p.Eta())< 2.4) {

      t_nEvent.eta1 = sortedGoodLeptons[PosLepton1].p.Eta();
      t_nEvent.pt1 = sortedGoodLeptons[PosLepton1].p.Pt();
      t_nEvent.phi1 = sortedGoodLeptons[PosLepton1].p.Phi();
      t_nEvent.ch1 = sortedGoodLeptons[PosLepton1].charge;
      t_nEvent.id1 = sortedGoodLeptons[PosLepton1].type; //??????
      for(int n = 0; n < nAnalysis; ++n) {
        t_nEvent.tagReco1[n] = sortedGoodLeptons[PosLepton1].tagReco[n]; //??????
        t_nEvent.probeReco1[n] = sortedGoodLeptons[PosLepton1].probeReco[n]; //??????
        t_nEvent.pprobeReco1[n] = sortedGoodLeptons[PosLepton1].pprobeReco[n]; //??????
        t_nEvent.tagIso1[n] = sortedGoodLeptons[PosLepton1].tagIso[n]; //??????
        t_nEvent.probeIso1[n] = sortedGoodLeptons[PosLepton1].probeIso[n]; //??????
        t_nEvent.pprobeIso1[n] = sortedGoodLeptons[PosLepton1].pprobeIso[n]; //??????
        t_nEvent.tagID1[n] = sortedGoodLeptons[PosLepton1].tagID[n]; //??????
        t_nEvent.probeID1[n] = sortedGoodLeptons[PosLepton1].probeID[n]; //??????
        t_nEvent.pprobeID1[n] = sortedGoodLeptons[PosLepton1].pprobeID[n]; //??????
        t_nEvent.tag1[n] = t_nEvent.tagReco1[n]*t_nEvent.tagIso1[n]*t_nEvent.tagID1[n]; 
        t_nEvent.probe1[n] = t_nEvent.probeReco1[n]*t_nEvent.probeIso1[n]*t_nEvent.probeID1[n];
        t_nEvent.pprobe1[n] = t_nEvent.pprobeReco1[n]*t_nEvent.pprobeIso1[n]*t_nEvent.pprobeID1[n];
      }
      t_nEvent.d0bs1 = sortedGoodLeptons[PosLepton1].d0BS;
      t_nEvent.dzbs1 = sortedGoodLeptons[PosLepton1].dzBS;
      t_nEvent.d0pv1 = sortedGoodLeptons[PosLepton1].d0PV;
      t_nEvent.dzpv1 = sortedGoodLeptons[PosLepton1].dzPV;
      if(sortedGoodLeptons[PosLepton1].pfindex != -1) {
        t_nEvent.pfeta1 = sortedGoodLeptons[PosLepton1].pfp.Eta();
        t_nEvent.pfpt1 = sortedGoodLeptons[PosLepton1].pfp.Pt();
        t_nEvent.pfphi1 = sortedGoodLeptons[PosLepton1].pfp.Phi();
      }
      t_nEvent.eta2 = sortedGoodLeptons[PosLepton2].p.Eta();
      t_nEvent.pt2 = sortedGoodLeptons[PosLepton2].p.Pt();
      t_nEvent.phi2 = sortedGoodLeptons[PosLepton2].p.Phi();
      t_nEvent.ch2 = sortedGoodLeptons[PosLepton2].charge;
      t_nEvent.id2 = sortedGoodLeptons[PosLepton2].type; //??????
      for(int n = 0; n < nAnalysis; ++n) {
        t_nEvent.tagReco2[n] = sortedGoodLeptons[PosLepton2].tagReco[n]; //??????
        t_nEvent.probeReco2[n] = sortedGoodLeptons[PosLepton2].probeReco[n]; //??????
        t_nEvent.pprobeReco2[n] = sortedGoodLeptons[PosLepton2].pprobeReco[n]; //??????
        t_nEvent.tagIso2[n] = sortedGoodLeptons[PosLepton2].tagIso[n]; //??????
        t_nEvent.probeIso2[n] = sortedGoodLeptons[PosLepton2].probeIso[n]; //??????
        t_nEvent.pprobeIso2[n] = sortedGoodLeptons[PosLepton2].pprobeIso[n]; //??????
        t_nEvent.tagID2[n] = sortedGoodLeptons[PosLepton2].tagID[n]; //??????
        t_nEvent.probeID2[n] = sortedGoodLeptons[PosLepton2].probeID[n]; //??????
        t_nEvent.pprobeID2[n] = sortedGoodLeptons[PosLepton2].pprobeID[n]; //??????
        t_nEvent.tag2[n] = t_nEvent.tagReco2[n]*t_nEvent.tagIso2[n]*t_nEvent.tagID2[n]; 
        t_nEvent.probe2[n] = t_nEvent.probeReco2[n]*t_nEvent.probeIso2[n]*t_nEvent.probeID2[n];
        t_nEvent.pprobe2[n] = t_nEvent.pprobeReco2[n]*t_nEvent.pprobeIso2[n]*t_nEvent.pprobeID2[n];
      }    
      t_nEvent.d0bs2 = sortedGoodLeptons[PosLepton2].d0BS;
      t_nEvent.dzbs2 = sortedGoodLeptons[PosLepton2].dzBS;
      t_nEvent.d0pv2 = sortedGoodLeptons[PosLepton2].d0PV;
      t_nEvent.dzpv2 = sortedGoodLeptons[PosLepton2].dzPV;
      if(sortedGoodLeptons[PosLepton2].pfindex != -1) {
        t_nEvent.pfeta2 = sortedGoodLeptons[PosLepton2].pfp.Eta();
        t_nEvent.pfpt2 = sortedGoodLeptons[PosLepton2].pfp.Pt();
        t_nEvent.pfphi2 = sortedGoodLeptons[PosLepton2].pfp.Phi();
      }
      t_nEvent.drl = sortedGoodLeptons[PosLepton1].p.DeltaR(sortedGoodLeptons[PosLepton2].p);
        

      negativePosition = PosLepton1;
      if(t_nEvent.ch2<0) negativePosition = PosLepton2;
      t_nEvent.etan = sortedGoodLeptons[negativePosition].p.Eta();
      t_nEvent.ptn = sortedGoodLeptons[negativePosition].p.Pt();
      t_nEvent.phin = sortedGoodLeptons[negativePosition].p.Phi();
      t_nEvent.chn = sortedGoodLeptons[negativePosition].charge;
      t_nEvent.idn = sortedGoodLeptons[negativePosition].type; //??????
      t_nEvent.d0bsn = sortedGoodLeptons[negativePosition].d0BS;
      t_nEvent.dzbsn = sortedGoodLeptons[negativePosition].dzBS;
      t_nEvent.d0pvn = sortedGoodLeptons[negativePosition].d0PV;
      t_nEvent.dzpvn = sortedGoodLeptons[negativePosition].dzPV;
      if(sortedGoodLeptons[negativePosition].pfindex != -1) {
        t_nEvent.pfetan = sortedGoodLeptons[negativePosition].pfp.Eta();
        t_nEvent.pfptn = sortedGoodLeptons[negativePosition].pfp.Pt();
        t_nEvent.pfphin = sortedGoodLeptons[negativePosition].pfp.Phi();
      }

      t_nEvent.mll=(sortedGoodLeptons[PosLepton2].p+sortedGoodLeptons[PosLepton1].p).M();
      t_nEvent.phi=(sortedGoodLeptons[PosLepton2].p+sortedGoodLeptons[PosLepton1].p).Phi();
      t_nEvent.pt=(sortedGoodLeptons[PosLepton2].p+sortedGoodLeptons[PosLepton1].p).Pt();
      if(sortedGoodLeptons[PosLepton2].pfindex != -1 && sortedGoodLeptons[PosLepton1].pfindex != -1) {
        t_nEvent.pfmll=(sortedGoodLeptons[PosLepton2].pfp+sortedGoodLeptons[PosLepton1].pfp).M();
        t_nEvent.pfphi=(sortedGoodLeptons[PosLepton2].pfp+sortedGoodLeptons[PosLepton1].pfp).Phi();
        t_nEvent.pfpt=(sortedGoodLeptons[PosLepton2].pfp+sortedGoodLeptons[PosLepton1].pfp).Pt();
        t_nEvent.pfdrl = sortedGoodLeptons[PosLepton1].pfp.DeltaR(sortedGoodLeptons[PosLepton2].pfp);
      }
      //t_nEvent.dphi=sortedGoodLeptons[PosLepton2].p.DeltaPhi(sortedGoodLeptons[PosLepton1].p);
      
      //t_nEvent.pfdphi=sortedGoodLeptons[PosLepton2].pfp.DeltaPhi(sortedGoodLeptons[PosLepton1].pfp);

    } else {
      
      if( isMC ) t_myTree->Fill();
      return;
      
    }

    
    float jetThreshold[nAnalysis] = {30, 40, 20};
    for(int n = 0; n < nAnalysis; n++) { 
      t_nEvent.pfJetGoodNum[n]=0;
      vector<t_lepton> pfGoodJets;
      float closest1=1000, closest2=1000, closestn=1000; //the dr distance between the leptons and jets
      for(int i =0 ; i<fTR->PF2PAT3NJets;i++) { // PF jet loop//killPF
	
        if(i==jMax){cout<<"max Num was reached"<<endl; return;}
	
	float jpt = fTR->PF2PAT3JPt[i];//killPF
	float jeta = fTR->PF2PAT3JEta[i];//killPF
	float jphi = fTR->PF2PAT3JPhi[i];//killPF
	float jpx = fTR->PF2PAT3JPx[i];//killPF
	float jpy = fTR->PF2PAT3JPy[i];//killPF
	float jpz = fTR->PF2PAT3JPz[i];//killPF
	float jenergy = fTR->PF2PAT3JE[i];//killPF
	//float jesC = fTR->PF2PAT3JEcorr[i];//killPF
	
	bool  isJetID = IsGoodBasicPFJet(i,false);
        if(n == 0) isJetID = IsGoodBasicPFJet(i, false); //JZB
	if(n == 1) isJetID = IsGoodBasicPFJet(i, false); //SS
	if(n == 2) isJetID = IsGoodBasicPFJetPAT3(i, 20.0, 3.0); //MT2
	
	TLorentzVector aJet(jpx,jpy,jpz,jenergy);
	TLorentzVector sumOfPFJets(0,0,0,0);
	// lepton-jet cleaning
	if ( fFullCleaning_ ) { 
	  // Remove jet close to any lepton
	  bool isClean(true);
	  for ( size_t ilep = 0; ilep<sortedGoodLeptons.size(); ++ilep )
	    if ( aJet.DeltaR(sortedGoodLeptons[ilep].p)<DRmax) isClean=false;
	  if ( !isClean ) continue;
	} else {
	  // Remove jet close to leptons from Z candidate
	  if(aJet.DeltaR(sortedGoodLeptons[PosLepton1].p)<DRmax)continue; 
	  if(aJet.DeltaR(sortedGoodLeptons[PosLepton2].p)<DRmax)continue;
	}
	
	// Keep jets over min. pt threshold
	if ( !(jpt>20) ) continue;  
	
	t_nEvent.pfHT[n]    += jpt;
	
	// Keep central jets
	if ( !(fabs(jeta)<2.6 ) ) continue;
	
	
	t_nEvent.pfGoodHT[n] += jpt;
	sumOfPFJets += aJet;
	
	t_lepton tmpLepton;
	tmpLepton.p = aJet;
	tmpLepton.charge = 0;
	tmpLepton.index = i;
	tmpLepton.type = -1;
	pfGoodJets.push_back(tmpLepton);
	
	if ( jpt > jetThreshold[n] ) {
	  t_nEvent.pfTightHT[n] += jpt;
	  if(isJetID>0) t_nEvent.pfJetGoodNumID[n]++;
	  t_nEvent.pfJetGoodNum[n]++;
	  float dr1 = aJet.DeltaR(sortedGoodLeptons[PosLepton1].p);
	  float dr2 = aJet.DeltaR(sortedGoodLeptons[PosLepton2].p);
	  float drn = aJet.DeltaR(sortedGoodLeptons[negativePosition].p);
	  if(dr1<closest1) closest1 = dr1;
	  if(dr2<closest2) closest2 = dr2;
	  if(drn<closestn) closestn = drn;
	}
	
      }
      
      t_nEvent.pfMET[n] = fTR->PFMET;
      t_nEvent.dr1[n] = closest1;
      t_nEvent.dr2[n] = closest2;
      t_nEvent.drn[n] = closestn;
    }
      t_myTree->Fill();	
  } else {
  
    if( isMC ) t_myTree->Fill();
   

  }

}


//__________________________________________________________________________
void RunEfficiency::End(){
  fHistFile->cd();	

  t_myTree->Write();

  fHistFile->Close();
}


//__________________________________________________________________________
template<class T>
std::string RunEfficiency::any2string(T i)
{
  std::ostringstream buffer;
  buffer << i;
  return buffer.str();
}



//__________________________________________________________________________
const bool RunEfficiency::IsCustomPfMu(const int index){
//VERY TEMPORARY !!!!
  // Basic muon cleaning and ID
  // Acceptance cuts
  if ( !(fTR->PfMu3Pt[index] > 5) )       return false;
  if ( !(fabs(fTR->PfMu3Eta[index])<2.4) ) return false;
  
/*  // Quality cuts
  if ( !fTR->fTR->PfMu3IsGMPT[index])        return false;
  if ( !fTR->MuIsGlobalMuon[index] )  return false;
  if ( !fTR->MuIsTrackerMuon[index] ) return false;
  
  // Hits
  if ( !(fTR->MuNTkHits[index] >= 11) )     return false;
  if ( !(fTR->MuNPxHits[index] > 0) )       return false;
  if ( !(fTR->MuNMatches[index] > 1) )      return false;
  
  // Vertex compatibility
  if ( !(fabs(fTR->MuD0PV[index]) < 0.02) ) return false;
  if ( !(fabs(fTR->MuDzPV[index]) < 1.0 ) ) return false;
  
  // Flat isolation below 20 GeV (only for synch.: we cut at 20...)
  double hybridIso = fTR->MuRelIso03[index]*fTR->MuPt[index]/std::max((float)20.,fTR->MuPt[index]);
  if ( !(hybridIso < 0.15) ) return false;
  
  if ( !(fTR->MuPtE[index]/fTR->MuPt[index] < 0.1) ) return false;
*/
  return true;
}

const bool RunEfficiency::IsCustomMu(const int index){

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


const bool RunEfficiency::IsCustomEl(const int index){

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


//__________________________________________________________________________
const bool RunEfficiency::IsCustomPfEl(const int index){

  // kinematic acceptance
  if(!(fTR->PfEl3Pt[index]>5) )return false;
  if(!(fabs(fTR->PfEl3Eta[index]) < 2.4) ) return false;
  if(!(fTR->PfElID95[index])) return false;
  //if(!(fTR->PfElID80[index])) return false;
/*
  if ( !(fTR->ElNumberOfMissingInnerHits[index] <= 1 ) ) return false;
  if ( !(fabs(fTR->ElD0PV[index]) < 0.04) ) return false;
  if ( !(fabs(fTR->ElDzPV[index]) < 1.0 ) ) return false;
  
  // Electron ID
  int elIDWP95 = fTR->ElIDsimpleWP95relIso[index];
  if (elIDWP95!=7) return false;
  
  // Flat isolation below 20 GeV (only for synch.)
  double hybridIso = fTR->ElRelIso03[index]
    *fTR->ElPt[index]/std::max((float)20.,fTR->ElPt[index]);
  if ( !(hybridIso < 0.15) ) return false;  
  
  //   // Other choices for electron ID
  //   if ( fTR->ElIDsimpleWP90relIso[index]!=7 ) return false;
  //   if ( fTR->ElIDsimpleWP80relIso[index]!=7 ) return false;
  //   if ( !(fTR->ElIDMva[index]>0.4) ) return false;
*/
  return true;
  
}


//_______________________________________________________________________
const int RunEfficiency::getPFMuIndex(const int recoIndex) {

  int muIndex = 0;
  float px= fTR->MuPx[recoIndex];
  float py= fTR->MuPy[recoIndex];
  float pz= fTR->MuPz[recoIndex];
  float energy =  fTR->MuE[recoIndex];
  TLorentzVector tmpVector(px,py,pz,energy);
  
  for(muIndex=0;muIndex<fTR->PfMu3NObjs;muIndex++) {
    float pfpx= fTR->PfMu3Px[muIndex];
    float pfpy= fTR->PfMu3Py[muIndex];
    float pfpz= fTR->PfMu3Pz[muIndex];
    float pfenergy =  fTR->PfMu3E[muIndex];
    //if(pfpx == 0 && pfpy == 0 && pfpz == 0) break;
    TLorentzVector tmpPfVector(pfpx,pfpy,pfpz,pfenergy);
    if(tmpPfVector.DeltaR(tmpVector)<0.05) return muIndex; //maybe to tight
  }
  return -1;
  
} 

//_______________________________________________________________________
const int RunEfficiency::getPFElIndex(const int recoIndex) {

  int elIndex = 0;
  float px= fTR->ElPx[recoIndex];
  float py= fTR->ElPy[recoIndex];
  float pz= fTR->ElPz[recoIndex];
  float energy =  fTR->ElE[recoIndex];
  TLorentzVector tmpVector(px,py,pz,energy);
  
  for(elIndex=0;elIndex<fTR->PfEl3NObjs;elIndex++) {
    float pfpx= fTR->PfEl3Px[elIndex];
    float pfpy= fTR->PfEl3Py[elIndex];
    float pfpz= fTR->PfEl3Pz[elIndex];
    float pfenergy =  fTR->PfEl3E[elIndex];
    //if(pfpx == 0 && pfpy == 0 && pfpz == 0) break;
    TLorentzVector tmpPfVector(pfpx,pfpy,pfpz,pfenergy);
    if(tmpPfVector.DeltaR(tmpVector)<0.05) return elIndex; //maybe to tight
  }
  return -1;

} 



const bool RunEfficiency::IsGoodBasicPFJetPAT3(int index, double ptcut, double absetacut){
        // Basic PF jet cleaning and ID cuts
        // cut at pt of ptcut (default = 30 GeV)
        // cut at abs(eta) of absetacut (default = 2.5)
        if(fTR->PF2PAT3JPt[index] < ptcut           ) return false;
        if(fabs(fTR->PF2PAT3JEta[index]) > absetacut) return false;
        if(fTR->PF2PAT3JIDLoose[index]    ==0       ) return false;
        if(fTR->PF2PAT3JScale[index]     < 0        ) return false;
        return true;
}


//__________________________________________________________________________
const bool RunEfficiency::IsCustomJet(const int index){
  // Basic Jet ID cuts (loose Jet ID)
  // See https://twiki.cern.ch/twiki/bin/view/CMS/JetID

  if ( !(fTR->CAJID_n90Hits[index] > 1) ) return false;
  if ( !(fTR->CAJID_HPD[index] < 0.98)  ) return false;
  
  if ( fabs(fTR->CAJEta[index])<2.6 ) {
    if ( !(fTR->CAJEMfrac[index] > 0.01)  ) return false;
  } else {
    if ( !(fTR->CAJEMfrac[index] > -0.9)  ) return false;
    if ( fTR->CAJPt[index] > 80 && !(fTR->CAJEMfrac[index]<1) ) return false;
  }
  
  return true;
}




//__________________________________________________________________________
const bool RunEfficiency::ElPassingTag(const int n, const int index){
  
  // kinematic acceptance
  if(!(fTR->ElPt[index]>5) )return false;
  if(!(fabs(fTR->ElEta[index]) < 2.4) ) return false;
  // Electron ID
  int elIDWP95 = fTR->ElIDsimpleWP95relIso[index];
  if (elIDWP95!=7) return false;

  //   // Other choices for electron ID
  //   if ( fTR->ElIDsimpleWP90relIso[index]!=7 ) return false;
  //   if ( fTR->ElIDsimpleWP80relIso[index]!=7 ) return false;
  //   if ( !(fTR->ElIDMva[index]>0.4) ) return false;
  return true;

}


//__________________________________________________________________________
const bool RunEfficiency::ElPassingRecoProbe(int n, int index){

   // kinematic acceptance
  if(!(fTR->ElPt[index]>5) )return false;
  if(!(fabs(fTR->ElEta[index]) < 2.4) ) return false; // Acceptance cuts
  
  if(n == 0) {
    if(!fTR->ElEcalDriven[index]) return false;
    if(fTR->ElCaloEnergy[index] < 10.) return false;
  } else if(n == 1) {
    if(!fTR->ElEcalDriven[index]) return false;
    if(fTR->ElCaloEnergy[index] < 10.) return false;
  } else if(n == 2) {
    if(!fTR->ElEcalDriven[index]) return false;
    if(fTR->ElCaloEnergy[index] < 10.) return false;
  }
  
  return true;

}

//__________________________________________________________________________
const bool RunEfficiency::ElPassingIDProbe(int n, int index){

  // Acceptance cuts
  if ( !(fTR->MuPt[index] > 5) )       return false;
  if ( !(fabs(fTR->MuEta[index])<2.4) ) return false;
  
  if(n == 0) {
   
  } else if(n == 1) {
 
  } else if(n == 2) {
    
  }

  return true;
}



//_________________________________________________________________________
const bool RunEfficiency::ElPassingTag(const int index) {

  // kinematic acceptance
  if(!(fTR->ElPt[index]>5) )return false;
  if(!(fabs(fTR->ElEta[index]) < 2.4) ) return false;
  // Electron ID
  int elIDWP95 = fTR->ElIDsimpleWP95relIso[index];
  if (elIDWP95!=7) return false;

  //   // Other choices for electron ID
  //   if ( fTR->ElIDsimpleWP90relIso[index]!=7 ) return false;
  //   if ( fTR->ElIDsimpleWP80relIso[index]!=7 ) return false;
  //   if ( !(fTR->ElIDMva[index]>0.4) ) return false;
  return true;
}

//____________________________________________________________________________
const bool RunEfficiency::ElPassingPProbe(const int index, int indexpf) {
  
  // kinematic acceptance
  if(!(fTR->ElPt[index]>5) )return false;
  if(!(fabs(fTR->ElEta[index]) < 2.4) ) return false;
  //if ( !(fTR->ElNumberOfMissingInnerHits[index] <= 1 ) ) return false;
  //if ( !(fabs(fTR->ElD0PV[index]) < 0.04) ) return false;
  //if ( !(fabs(fTR->ElDzPV[index]) < 1.0 ) ) return false;
  if ( indexpf == -1) return false;
  /*
  // Electron ID
  int elIDWP95 = fTR->ElIDsimpleWP95relIso[index];
  if (elIDWP95!=7) return false;

  // Flat isolation below 20 GeV (only for synch.)
  double hybridIso = fTR->ElRelIso03[index]
    *fTR->ElPt[index]/std::max((float)20.,fTR->ElPt[index]);
  if ( !(hybridIso < 0.15) ) return false;

  //   // Other choices for electron ID
  //   if ( fTR->ElIDsimpleWP90relIso[index]!=7 ) return false;
  //   counters[EL].fill("... passes WP90 ID");
  //   if ( fTR->ElIDsimpleWP80relIso[index]!=7 ) return false;
  //   counters[EL].fill("... passes WP80 ID");
  //   if ( !(fTR->ElIDMva[index]>0.4) ) return false;
  //   counters[EL].fill("... MVA>0.4");
  */
  if( ! IsCustomPfEl(indexpf) ) return false;
  

  return true;
}

//__________________________________________________________________________
bool RunEfficiency::GeneratorInfo(void) {
 

 // Try to find an Z->ll pair inside the acceptance
  double minPt = 5.;
  double mllCut = 30.;
  double maxEta = 2.4;
  double minJPt = 20;
  double maxJEta = 3.0;
  bool saveMC = false;
  // First, look for leptons in acceptance
  vector<t_lepton> gLeptons;
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
        t_lepton tmpLepton;
        tmpLepton.p      = tmpVector;
        tmpLepton.charge = fTR->GenLeptonID[gIndex]/abs(fTR->GenLeptonID[gIndex]);
        tmpLepton.index  = gIndex;
        tmpLepton.type   = fTR->GenLeptonID[gIndex];
        tmpLepton.genPt  = tmpVector.Pt();
        if ( fTR->GenLeptonMID[gIndex] ==23 ) gLeptons.push_back(tmpLepton); // WW study
      }         
  }
  // Gen leptons are not sorted by Pt...
  vector<t_lepton> sortedGLeptons = sortLeptonsByPt(gLeptons);
  
  // Store actual number of leptons passing selection
  t_nEvent.genNleptons = gLeptons.size();
  
  // Now fill information
  TLorentzVector GenMETvector(fTR->GenMETpx,fTR->GenMETpy,0,0);
  t_nEvent.genMET     = fTR->GenMET;
  
  // Number of good jets
  t_nEvent.genNjets = 0;
  for ( int jIndex=0; jIndex<fTR->NGenJets; ++jIndex) {
    if ( fTR->GenJetPt[jIndex]<minJPt ) continue;
    if ( fabs(fTR->GenJetEta[jIndex])>maxJEta ) continue;
    ++t_nEvent.genNjets;
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
    //if(i1 == sortedGLeptons.size()-1) return saveMC;
  }
 

  
 
  
  float dr1 = 1000, dr2 = 1000, drN = 1000;
  
  if(sortedGLeptons.size()>0) {
    t_nEvent.genPt1     = sortedGLeptons[i1].p.Pt();
    t_nEvent.genId1     = sortedGLeptons[i1].type;
    t_nEvent.genEta1    = sortedGLeptons[i1].p.Eta();
    t_nEvent.genMID1     = fTR->GenLeptonMID[sortedGLeptons[i1].index];
    if(sortedGLeptons.size()>1) {
      // Number of good jets
      /*for ( int jIndex=0; jIndex<fTR->NGenJets; ++jIndex) {
        if ( fTR->GenJetPt[jIndex]<minJPt ) continue;
        if ( fabs(fTR->GenJetEta[jIndex])>maxJEta ) continue;
        TLorentzVector theJet(fTR->GenJetPx[jIndex], fTR->GenJetPy[jIndex], fTR->GenJetPz[jIndex], fTR->GenJetE[jIndex]);
        float dr1_, dr2_;
        dr1_ = theJet.DeltaR(sortedGLeptons[i1].p);
        dr2_ = theJet.DeltaR(sortedGLeptons[i2].p);
        if(dr1_ < dr1) dr1 = dr1_;  
        if(dr2_ < dr2) dr2 = dr2_;  
      }
      */
      TLorentzVector genZvector = sortedGLeptons[i1].p + sortedGLeptons[i2].p;
      t_nEvent.genPt2     = sortedGLeptons[i2].p.Pt();
      t_nEvent.genId2     = sortedGLeptons[i2].type;
      t_nEvent.genEta2    = sortedGLeptons[i2].p.Eta();
      t_nEvent.genZPt     = genZvector.Pt();
      t_nEvent.genMll     = genZvector.M();
      t_nEvent.genMID2     = fTR->GenLeptonMID[sortedGLeptons[i2].index];
      saveMC = true;
      if(t_nEvent.genId2<0) {
        t_nEvent.genPtN     = t_nEvent.genPt2;
        t_nEvent.genEtaN    = t_nEvent.genEta2;
        t_nEvent.genDRN     = dr2;  
      } else {
        t_nEvent.genPtN     = t_nEvent.genPt1;
        t_nEvent.genEtaN    = t_nEvent.genEta1; 
        t_nEvent.genDRN     = dr1;  
      }
    }
  }
  
  return saveMC;
}


//********************************* FUNCTIONS THAT ARE STILL MISSING !!!!

const bool RunEfficiency::ElPassingIsoTag(int a, int b) {
cout << "NOT DEFINED!" << endl;
}

const bool RunEfficiency::ElPassingIsoProbe(int a, int b) {
cout << "NOT DEFINED!" << endl;
}

const bool RunEfficiency::ElPassingIDTag(int a, int b) {
cout << "NOT DEFINED!" << endl;
}

const bool RunEfficiency::ElPassingIDPProbe(int a,int b,int c) {
cout << "NOT DEFINED!" << endl;
}

const bool RunEfficiency::ElPassingRecoTag(int a, int b) {
cout << "NOT DEFINED!" << endl;
}
const bool RunEfficiency::ElPassingRecoPProbe(int a, int b) {
cout << "NOT DEFINED!" << endl;
}
const bool RunEfficiency::ElPassingIsoPProbe(int a, int b) {
cout << "NOT DEFINED!" << endl;
}


//************************************************************************************

/*

Steps: 		Tag	Probe	Passing Probe
(basic)		1	
Reco		-	2	3
ID		-	3	4
Isolation	-	4	5


1) Tag
2) Reco Probe
3) Passing Reco Probe
- 3) ID Probe
4) Passing ID Probe
- 4) Iso Probe
5) Iso Passing Probe

*/
 
///
/// 1
///

const bool RunEfficiency::MuPassingTag(int n, int index){

  // Acceptance cuts
  if(n!=2) {
	if ( !(fTR->MuPt[index] > 5) )       return false;
	if ( !(fabs(fTR->MuEta[index])<2.4) ) return false;
  } else {
	if(getPFMuIndex(index)==-1) return false;
  }
  /*
  if(n == 0) {
    // Quality cuts
    if ( !fTR->MuIsGMPT[index] )        return false;
    if ( !fTR->MuIsGlobalMuon[index] )  return false;
    if ( !fTR->MuIsTrackerMuon[index] ) return false;
  } else if(n == 1) {
    if ( !fTR->MuIsGMPT[index] )        return false;
    if ( !fTR->MuIsGlobalMuon[index] )  return false;
    if ( !fTR->MuIsTrackerMuon[index] ) return false;
  } else if(n == 2) {
    if ( !fTR->MuIsGMPT[index] )        return false;
    if ( !fTR->MuIsGlobalMuon[index] )  return false;
    if ( !fTR->MuIsTrackerMuon[index] ) return false;
  }
*/
  return true;
}

const bool RunEfficiency::MuPassingRecoTag(int n, int index){
  return RunEfficiency::MuPassingTag(n,index);
}

const bool RunEfficiency::MuPassingIDTag(int n, int index) {
  return RunEfficiency::MuPassingTag(n,index);
}

const bool RunEfficiency::MuPassingIsoTag(int n, int index) {
  return RunEfficiency::MuPassingTag(n,index);
}

///--------------------------------------------------------------------------------
///
/// 2
///

const bool RunEfficiency::MuPassingRecoProbe(int n, int index){

  // Acceptance cuts
  if ( !(fTR->MuPt[index] > 5) )       return false;
  if ( !(fabs(fTR->MuEta[index])<2.4) ) return false;
  

  return true;
}

///--------------------------------------------------------------------------------
///
/// 3
///

const bool RunEfficiency::MuPassingRecoPProbe(const int n, const int index, const int pfindex){
  
  // Acceptance cuts
  if ( !(fTR->MuPt[index] > 5) )       return false;
  if ( !(fabs(fTR->MuEta[index])<2.4) ) return false;

  if(n == 0) {//JZB
	if ( !(fTR->MuNTkHits[index] >= 11) )     return false;
	if ( !(fTR->MuNPxHits[index] > 0) )       return false;
	if ( !(fTR->MuNMatches[index] > 1) )      return false;
	if ( !(fabs(fTR->MuD0PV[index]) < 0.02) ) return false;
	if ( !(fabs(fTR->MuDzPV[index]) < 1.0 ) ) return false;
  } else if(n == 1) {//SameSign
	if(fTR->MuNChi2[index] > 10)   return false;
	if(fTR->MuPtE[index]/fTR->MuPt[index] > 0.1) return false;
	if(fTR->MuNTkHits[index] < 11) return false;
	if(fTR->MuNMuHits[index] < 1)  return false;
	if(fabs(fTR->MuD0PV[index]) > 0.02)    return false;
	if(fabs(fTR->MuDzPV[index]) > 1.00)    return false;
	if(fTR->MuIso03EMVetoEt[index] > 4.0)  return false;
	if(fTR->MuIso03HadVetoEt[index] > 6.0) return false;
  } else if(n == 2) {//MT2
//	if ( !fTR->MuIsTrackerMuon[index] ) return false;
	if(indexpf == -1) return false;
	cout << "don't know exactly what to do for mt2!" << endl;
  }
  return true;
}

/// 
/// 3a
///

const bool RunEfficiency::MuPassingIDProbe(int n, int index){
  return RunEfficiency::MuPassingRecoPProbe(n,index,getPFMuIndex(index));
} 

///--------------------------------------------------------------------------------
///
/// 4
///


const bool RunEfficiency::MuPassingIDPProbe(int n, int index, int pfindex){
  if(!RunEfficiency::MuPassingIDProbe(n,index)) return false;
  if(n==0) // JZB
  {
	if ( !fTR->MuIsGMPT[index] )        return false;
	if ( !fTR->MuIsGlobalMuon[index] )  return false;
	if ( !fTR->MuIsTrackerMuon[index] ) return false;
  } else if(n==1) // same sign
  {
	if(fTR->MuIsGlobalMuon[index] == 0)  return false;
	if(fTR->MuIsTrackerMuon[index] == 0) return false;
  } else if(n==1) // mt2
  {
	if(indexpf == -1) return false;
	cout << "don't know exactly what to do for mt2!" << endl;
  }
  return true;
}

///
/// 4a
///

const bool RunEfficiency::MuPassingIsoProbe(int n, int index){
  return RunEfficiency::MuPassingIDPProbe(n, index, getPFMuIndex(index));
} 

///--------------------------------------------------------------------------------
///
/// 5
///

const bool RunEfficiency::MuPassingIsoPProbe(int n, int index, int pfindex){
  RunEfficiency::MuPassingIDPProbe(n,index,pfindex);
  if(n==0) // JZB
  {
	double hybridIso = fTR->MuRelIso03[index]*fTR->MuPt[index]/std::max((float)20.,fTR->MuPt[index]);
	if ( !(hybridIso < 0.15) ) return false;
	if ( !(fTR->MuPtE[index]/fTR->MuPt[index] < 0.1) ) return false;
  } else if(n==1) // same sign
  {
	if(fTR->MuRelIso03[index] > 0.15) return false;
  } else if(n==1) // same sign
  {
	if(indexpf == -1) return false;
	cout << "don't know exactly what to do for mt2!" << endl;
  }
  return true;
} 

/*

Steps: 		Tag	Probe	Passing Probe
(basic)		1	
Reco		-	2	3
ID		-	3	4
Isolation	-	4	5


1) Tag
2) Reco Probe
3) Passing Reco Probe
- 3) ID Probe
4) Passing ID Probe
- 4) Iso Probe
5) Iso Passing Probe

*/