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
#include "TagNProbeDefsEle.cc"
#include "TagNProbeDefsMu.cc"


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
  float dr1jzb,dr1ss,dr1mt2;
  float dr2jzb,dr2ss,dr2mt2;
  float drnjzb,drnss,drnmt2;
  int ch1;
  int ch2;
  int chn;
  int pfch1;
  int pfch2;
  int pfchn;
  int id1;
  int id2;
  int idn;
  
  int probeReco1jzb,probeReco1ss,probeReco1mt2;
  int pprobeReco1jzb,pprobeReco1ss,pprobeReco1mt2;
  int probeReco2jzb,probeReco2ss,probeReco2mt2;
  int pprobeReco2jzb,pprobeReco2ss,pprobeReco2mt2;
  
  int probeIso1jzb,probeIso1ss,probeIso1mt2;
  int pprobeIso1jzb,pprobeIso1ss,pprobeIso1mt2;
  int probeIso2jzb,probeIso2ss,probeIso2mt2;
  int pprobeIso2jzb,pprobeIso2ss,pprobeIso2mt2;
  
  int probeID1jzb,probeID1ss,probeID1mt2;
  int pprobeID1jzb,pprobeID1ss,pprobeID1mt2;
  int probeID2jzb,probeID2ss,probeID2mt2;
  int pprobeID2jzb,pprobeID2ss,pprobeID2mt2;

  int tag1jzb,tag1ss,tag1mt2;
  int probe1jzb,probe1ss,probe1mt2;
  int pprobe1jzb,pprobe1ss,pprobe1mt2;
  int tag2jzb,tag2ss,tag2mt2;
  int probe2jzb,probe2ss,probe2mt2;
  int pprobe2jzb,pprobe2ss,pprobe2mt2;

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
  int pfJetGoodNumjzb,pfJetGoodNumss,pfJetGoodNummt2;
  int pfJetGoodNumIDjzb,pfJetGoodNumIDss,pfJetGoodNumIDmt2;
  float pfMETjzb,pfMETss,pfMETmt2;
  float pfHTjzb,pfHTss,pfHTmt2;
  float pfGoodHTjzb,pfGoodHTss,pfGoodHTmt2;
  float pfTightHTjzb,pfTightHTss,pfTightHTmt2;
  
  int eventNum;
  int runNum;
  int lumi;
  int goodVtx;
  int numVtx;
  
  int passedee_triggersjzb,passedee_triggersss,passedee_triggersmt2;
  int passedmm_triggersjzb,passedmm_triggersss,passedmm_triggersmt2;
  int passedem_triggersjzb,passedem_triggersss,passedem_triggersmt2;
 
  
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

    tag1jzb = 0;tag1ss = 0;tag1mt2 = 0; 
    probe1jzb = 0;probe1ss = 0;probe1mt2 = 0; 
    pprobe1jzb = 0;pprobe1ss = 0;pprobe1mt2 = 0;
    tag2jzb = 0; tag2ss = 0; tag2mt2 = 0; 
    probe2jzb = 0; probe2ss = 0; probe2mt2 = 0; 
    pprobe2jzb = 0;pprobe2ss = 0;pprobe2mt2 = 0;
  
    probeReco1jzb = 0; probeReco1ss = 0; probeReco1mt2 = 0; 
    pprobeReco1jzb = 0;pprobeReco1ss = 0;pprobeReco1mt2 = 0;
    probeReco2jzb = 0; probeReco2ss = 0; probeReco2mt2 = 0; 
    pprobeReco2jzb = 0;pprobeReco2ss = 0;pprobeReco2mt2 = 0;
  
    probeIso1jzb = 0; probeIso1ss = 0; probeIso1mt2 = 0; 
    pprobeIso1jzb = 0;pprobeIso1ss = 0;pprobeIso1mt2 = 0;
    probeIso2jzb = 0; probeIso2ss = 0; probeIso2mt2 = 0; 
    pprobeIso2jzb = 0;pprobeIso2ss = 0;pprobeIso2mt2 = 0;
    
    probeID1jzb = 0; probeID1ss = 0; probeID1mt2 = 0; 
    pprobeID1jzb = 0;pprobeID1ss = 0;pprobeID1mt2 = 0;
    probeID2jzb = 0; probeID2ss = 0; probeID2mt2 = 0; 
    pprobeID2jzb = 0;pprobeID2ss = 0;pprobeID2mt2 = 0;

    dr1jzb = 0; dr1ss = 0; dr1mt2 = 0; 
    dr2jzb = 0; dr2ss = 0; dr2mt2 = 0; 
    drnjzb = 0; drnss = 0; drnmt2 = 0; 
  
    pfJetGoodNumjzb = 0; pfJetGoodNumss = 0; pfJetGoodNummt2 = 0; 
    pfJetGoodNumIDjzb = 0;pfJetGoodNumIDss = 0;pfJetGoodNumIDmt2 = 0;
    pfHTjzb = 0; pfHTss = 0; pfHTmt2 = 0; 
    pfGoodHTjzb = 0; pfGoodHTss = 0; pfGoodHTmt2 = 0; 
    pfTightHTjzb = 0;pfTightHTss = 0;pfTightHTmt2 = 0;

    passedee_triggersjzb = 0; passedee_triggersss = 0; passedee_triggersmt2 = 0; 
    passedmm_triggersjzb = 0; passedmm_triggersss = 0;passedmm_triggersmt2 = 0;  
    passedem_triggersjzb = 0;passedem_triggersss = 0;passedem_triggersmt2 = 0;

    pfMETjzb=0;pfMETss=0;pfMETmt2=0;

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
  t_myTree->Branch("dr1jzb",&t_nEvent.dr1jzb,"dr1jzb/F");
  t_myTree->Branch("dr1ss",&t_nEvent.dr1ss,"dr1ss/F");
  t_myTree->Branch("dr1mt2",&t_nEvent.dr1mt2,"dr1mt2/F");
  t_myTree->Branch("dr2jzb",&t_nEvent.dr2jzb,"dr2jzb/F");
  t_myTree->Branch("dr2ss",&t_nEvent.dr2ss,"dr2ss/F");
  t_myTree->Branch("dr2mt2",&t_nEvent.dr2mt2,"dr2mt2/F");
  t_myTree->Branch("drnjzb",&t_nEvent.drnjzb,"drnjzb/F");
  t_myTree->Branch("drnss",&t_nEvent.drnss,"drnss/F");
  t_myTree->Branch("drnmt2",&t_nEvent.drnmt2,"drnmt2/F");
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
  t_myTree->Branch("pfMETjzb",&t_nEvent.pfMETjzb,"pfMETjzb/F");
  t_myTree->Branch("pfMETss",&t_nEvent.pfMETss,"pfMETss/F");
  t_myTree->Branch("pfMETmt2",&t_nEvent.pfMETmt2,"pfMETmt2/F");
  t_myTree->Branch("pfHTjzb",&t_nEvent.pfHTjzb,"pfHTjzb/F");
  t_myTree->Branch("pfHTss",&t_nEvent.pfHTss,"pfHTss/F");
  t_myTree->Branch("pfHTmt2",&t_nEvent.pfHTmt2,"pfHTmt2/F");
  t_myTree->Branch("pfGoodHTjzb", &t_nEvent.pfGoodHTjzb,"pfGoodHTjzb/F");
  t_myTree->Branch("pfGoodHTss", &t_nEvent.pfGoodHTss,"pfGoodHTss/F");
  t_myTree->Branch("pfGoodHTmt2", &t_nEvent.pfGoodHTmt2,"pfGoodHTmt2/F");
  t_myTree->Branch("pfTightHTjzb", &t_nEvent.pfTightHTjzb,"pfTightHTjzb/F");
  t_myTree->Branch("pfTightHTss", &t_nEvent.pfTightHTss,"pfTightHTss/F");
  t_myTree->Branch("pfTightHTmt2", &t_nEvent.pfTightHTmt2,"pfTightHTmt2/F");
  t_myTree->Branch("eventNum",&t_nEvent.eventNum,"eventNum/I");
  t_myTree->Branch("runNum",&t_nEvent.runNum,"runNum/I");
  t_myTree->Branch("lumi",&t_nEvent.lumi,"lumi/I");
  t_myTree->Branch("goodVtx",&t_nEvent.goodVtx,"goodVtx/I");
  t_myTree->Branch("numVtx",&t_nEvent.numVtx,"numVtx/I");
  t_myTree->Branch("pfJetGoodNumjzb",&t_nEvent.pfJetGoodNumjzb,"pfJetGoodNumjzb/I");
  t_myTree->Branch("pfJetGoodNumss",&t_nEvent.pfJetGoodNumss,"pfJetGoodNumss/I");
  t_myTree->Branch("pfJetGoodNummt2",&t_nEvent.pfJetGoodNummt2,"pfJetGoodNummt2/I");
  t_myTree->Branch("pfJetGoodNumIDjzb",&t_nEvent.pfJetGoodNumIDjzb,"pfJetGoodNumIDjzb/I");
  t_myTree->Branch("pfJetGoodNumIDss",&t_nEvent.pfJetGoodNumIDss,"pfJetGoodNumIDss/I");
  t_myTree->Branch("pfJetGoodNumIDmt2",&t_nEvent.pfJetGoodNumIDmt2,"pfJetGoodNumIDmt2/I");
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
  t_myTree->Branch("passedee_triggersjzb", &t_nEvent.passedee_triggersjzb,"passedee_triggersjzb/I");
  t_myTree->Branch("passedee_triggersss",  &t_nEvent.passedee_triggersss,"passedee_triggersss/I");
  t_myTree->Branch("passedee_triggersmt2", &t_nEvent.passedee_triggersmt2,"passedee_triggersmt2/I");
  t_myTree->Branch("passedmm_triggersjzb", &t_nEvent.passedmm_triggersjzb,"passedmm_triggersjzb/I");
  t_myTree->Branch("passedmm_triggersss",  &t_nEvent.passedmm_triggersss,"passedmm_triggersss/I");
  t_myTree->Branch("passedmm_triggersmt2", &t_nEvent.passedmm_triggersmt2,"passedmm_triggersmt2/I");
  t_myTree->Branch("passedem_triggersjzb", &t_nEvent.passedem_triggersjzb,"passedem_triggersjzb/I");
  t_myTree->Branch("passedem_triggersss",  &t_nEvent.passedem_triggersss,"passedem_triggersss/I");
  t_myTree->Branch("passedem_triggersmt2", &t_nEvent.passedem_triggersmt2,"passedem_triggersmt2/I");
  t_myTree->Branch("tag1jzb", &t_nEvent.tag1jzb,"tag1jzb/I");
  t_myTree->Branch("tag1ss", &t_nEvent.tag1ss,"tag1ss/I");
  t_myTree->Branch("tag1mt2", &t_nEvent.tag1mt2,"tag1mt2/I");
  t_myTree->Branch("probe1jzb", &t_nEvent.probe1jzb,"probe1jzb/I");
  t_myTree->Branch("probe1ss", &t_nEvent.probe1ss,"probe1ss/I");
  t_myTree->Branch("probe1mt2", &t_nEvent.probe1mt2,"probe1mt2/I");
  t_myTree->Branch("pprobe1jzb", &t_nEvent.pprobe1jzb,"pprobe1jzb/I");
  t_myTree->Branch("pprobe1ss", &t_nEvent.pprobe1ss,"pprobe1ss/I");
  t_myTree->Branch("pprobe1mt2", &t_nEvent.pprobe1mt2,"pprobe1mt2/I");
  t_myTree->Branch("tag2jzb", &t_nEvent.tag2jzb,"tag2jzb/I");
  t_myTree->Branch("tag2ss", &t_nEvent.tag2ss,"tag2ss/I");
  t_myTree->Branch("tag2mt2", &t_nEvent.tag2mt2,"tag2mt2/I");
  t_myTree->Branch("probe2jzb", &t_nEvent.probe2jzb,"probe2jzb/I");
  t_myTree->Branch("probe2ss", &t_nEvent.probe2ss,"probe2ss/I");
  t_myTree->Branch("probe2mt2", &t_nEvent.probe2mt2,"probe2mt2/I");
  t_myTree->Branch("pprobe2jzb", &t_nEvent.pprobe2jzb,"pprobe2jzb/I");
  t_myTree->Branch("pprobe2ss", &t_nEvent.pprobe2ss,"pprobe2ss/I");
  t_myTree->Branch("pprobe2mt2", &t_nEvent.pprobe2mt2,"pprobe2mt2/I");
  t_myTree->Branch("probeReco1jzb", &t_nEvent.probeReco1jzb,"probeReco1jzb/I");
  t_myTree->Branch("probeReco1ss", &t_nEvent.probeReco1ss,"probeReco1ss/I");
  t_myTree->Branch("probeReco1mt2", &t_nEvent.probeReco1mt2,"probeReco1mt2/I");
  t_myTree->Branch("probeReco2jzb", &t_nEvent.probeReco2jzb,"probeReco2jzb/I");
  t_myTree->Branch("probeReco2ss", &t_nEvent.probeReco2ss,"probeReco2ss/I");
  t_myTree->Branch("probeReco2mt2", &t_nEvent.probeReco2mt2,"probeReco2mt2/I");
  t_myTree->Branch("pprobeReco2jzb", &t_nEvent.pprobeReco2jzb,"pprobeReco2jzb/I");
  t_myTree->Branch("pprobeReco2ss", &t_nEvent.pprobeReco2ss,"pprobeReco2ss/I");
  t_myTree->Branch("pprobeReco2mt2", &t_nEvent.pprobeReco2mt2,"pprobeReco2mt2/I");
  t_myTree->Branch("probeIso1jzb", &t_nEvent.probeIso1jzb,"probeIso1jzb/I");
  t_myTree->Branch("probeIso1ss", &t_nEvent.probeIso1ss,"probeIso1ss/I");
  t_myTree->Branch("probeIso1mt2", &t_nEvent.probeIso1mt2,"probeIso1mt2/I");
  t_myTree->Branch("pprobeIso1jzb", &t_nEvent.pprobeIso1jzb,"pprobeIso1jzb/I");
  t_myTree->Branch("pprobeIso1ss", &t_nEvent.pprobeIso1ss,"pprobeIso1ss/I");
  t_myTree->Branch("pprobeIso1mt2", &t_nEvent.pprobeIso1mt2,"pprobeIso1mt2/I");
  t_myTree->Branch("probeIso2jzb", &t_nEvent.probeIso2jzb,"probeIso2jzb/I");
  t_myTree->Branch("probeIso2ss", &t_nEvent.probeIso2ss,"probeIso2ss/I");
  t_myTree->Branch("probeIso2mt2", &t_nEvent.probeIso2mt2,"probeIso2mt2/I");
  t_myTree->Branch("pprobeIso2jzb", &t_nEvent.pprobeIso2jzb,"pprobeIso2jzb/I");
  t_myTree->Branch("pprobeIso2ss", &t_nEvent.pprobeIso2ss,"pprobeIso2ss/I");
  t_myTree->Branch("pprobeIso2mt2", &t_nEvent.pprobeIso2mt2,"pprobeIso2mt2/I");
  t_myTree->Branch("probeID1jzb", &t_nEvent.probeID1jzb,"probeID1jzb/I");
  t_myTree->Branch("probeID1ss", &t_nEvent.probeID1ss,"probeID1ss/I");
  t_myTree->Branch("probeID1mt2", &t_nEvent.probeID1mt2,"probeID1mt2/I");
  t_myTree->Branch("pprobeID1jzb", &t_nEvent.pprobeID1jzb,"pprobeID1jzb/I");
  t_myTree->Branch("pprobeID1ss", &t_nEvent.pprobeID1ss,"pprobeID1ss/I");
  t_myTree->Branch("pprobeID1mt2", &t_nEvent.pprobeID1mt2,"pprobeID1mt2/I");
  t_myTree->Branch("probeID2jzb", &t_nEvent.probeID2jzb,"probeID2jzb/I");
  t_myTree->Branch("probeID2ss", &t_nEvent.probeID2ss,"probeID2ss/I");
  t_myTree->Branch("probeID2mt2", &t_nEvent.probeID2mt2,"probeID2mt2/I");
  t_myTree->Branch("pprobeID2jzb", &t_nEvent.pprobeID2jzb,"pprobeID2jzb/I");
  t_myTree->Branch("pprobeID2ss", &t_nEvent.pprobeID2ss,"pprobeID2ss/I");
  t_myTree->Branch("pprobeID2mt2", &t_nEvent.pprobeID2mt2,"pprobeID2mt2/I");
  
}



//________________________________________________________________________________
vector<t_lepton> RunEfficiency::sortLeptonsByPt(vector<t_lepton>& leptons) {
  
  vector<t_lepton> theLep = leptons;
  sort (theLep.begin(), theLep.end(), momentumComparator);
  return theLep;  
  
}

const bool RunEfficiency::passAnyMT2Trigger() {
// HT
if ( GetHLTResult("HLT_HT150_v3") )        return true;
if ( GetHLTResult("HLT_HT160_v2") )        return true;
if ( GetHLTResult("HLT_HT200_v2") )        return true;
if ( GetHLTResult("HLT_HT200_v3") )        return true;
if ( GetHLTResult("HLT_HT240_v2") )        return true;
if ( GetHLTResult("HLT_HT250_v2") )        return true;
if ( GetHLTResult("HLT_HT250_v3") )        return true;
if ( GetHLTResult("HLT_HT260_v2") )        return true;
if ( GetHLTResult("HLT_HT300_v2") )        return true;
if ( GetHLTResult("HLT_HT300_v3") )        return true;
if ( GetHLTResult("HLT_HT300_v4") )        return true;
if ( GetHLTResult("HLT_HT300_v5") )        return true;
if ( GetHLTResult("HLT_HT350_v2") )        return true;
if ( GetHLTResult("HLT_HT350_v3") )        return true;
if ( GetHLTResult("HLT_HT350_v4") )        return true;
if ( GetHLTResult("HLT_HT360_v2") )        return true;
if ( GetHLTResult("HLT_HT400_v2") )        return true;
if ( GetHLTResult("HLT_HT400_v3") )        return true;
if ( GetHLTResult("HLT_HT400_v4") )        return true;
if ( GetHLTResult("HLT_HT440_v2") )        return true;
if ( GetHLTResult("HLT_HT450_v2") )        return true;
if ( GetHLTResult("HLT_HT450_v3") )        return true;
if ( GetHLTResult("HLT_HT450_v4") )        return true;
if ( GetHLTResult("HLT_HT500_v2") )        return true;
if ( GetHLTResult("HLT_HT500_v3") )        return true;
if ( GetHLTResult("HLT_HT500_v4") )        return true;
if ( GetHLTResult("HLT_HT550_v2") )        return true;
if ( GetHLTResult("HLT_HT550_v3") )        return true;
if ( GetHLTResult("HLT_HT550_v4") )        return true;
if ( GetHLTResult("HLT_HT550_v5") )        return true;
// MHT_HT
if ( GetHLTResult("HLT_HT250_MHT60_v2") )        return true;
if ( GetHLTResult("HLT_HT250_MHT60_v3") )        return true;
if ( GetHLTResult("HLT_HT250_MHT60_v4") )        return true;
if ( GetHLTResult("HLT_HT250_MHT70_v1") )        return true;
if ( GetHLTResult("HLT_HT260_MHT60_v2") )        return true;
if ( GetHLTResult("HLT_HT300_MHT75_v4") )        return true;
if ( GetHLTResult("HLT_HT300_MHT75_v5") )        return true;
// QuadJet
if ( GetHLTResult("HLT_QuadJet50_BTagIP_v1") )        return true;
if ( GetHLTResult("HLT_QuadJet50_Jet40_v1") )        return true;
// Muons
if ( GetHLTResult("HLT_DoubleMu3_HT160_v2") )        return true;
if ( GetHLTResult("HLT_DoubleMu3_v3") )        return true;
if ( GetHLTResult("HLT_Mu8_Jet40_v2") )        return true;
return false;
}



//________________________________________________________________________________
const bool RunEfficiency::passElTriggers(int iAnalysis) {
  if(iAnalysis==0) {//JZB
	if ( GetHLTResult("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v1") )        return true;
	if ( GetHLTResult("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v2") )        return true;
	if ( GetHLTResult("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v3") )        return true;
	if ( GetHLTResult("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v4") )        return true;
	if ( GetHLTResult("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v5") )        return true;
  }
  if(iAnalysis==1) {//same sign
	if ( GetHLTResult("HLT_ELE8_JET40") )        				return true;
	if ( GetHLTResult("HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v1") )        	return true;
	if ( GetHLTResult("HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v2") )        	return true;
	if ( GetHLTResult("HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v3") )        	return true;
	if ( GetHLTResult("HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v4") )        	return true;
	if ( GetHLTResult("HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v5") )        	return true;
	if ( GetHLTResult("HLT_ELE17_ELE8") )        				return true;
	if ( GetHLTResult("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v1") )        return true;
	if ( GetHLTResult("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v2") )        return true;
	if ( GetHLTResult("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v3") )        return true;
	if ( GetHLTResult("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v4") )        return true;
	if ( GetHLTResult("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v5") )        return true;
	if ( GetHLTResult("HLT_ELE17_ELE8_TIGHT") )        			return true;
	if ( GetHLTResult("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v1") )        return true;
	if ( GetHLTResult("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v2") )        return true;
	if ( GetHLTResult("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v3") )        return true;
	if ( GetHLTResult("HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v1") )        return true;
	if ( GetHLTResult("HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v2") )        return true;
	if ( GetHLTResult("HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v3") )        return true;
	if ( GetHLTResult("HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v4") )        return true;
	if ( GetHLTResult("HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v5") )        return true;
	if ( GetHLTResult("HLT_DOUBLEELE8_HT160") )        			return true;
	if ( GetHLTResult("HLT_DoubleEle8_CaloIdL_TrkIdVL_HT160_v1") )        	return true;
	if ( GetHLTResult("HLT_DoubleEle8_CaloIdL_TrkIdVL_HT160_v2") )        	return true;
	if ( GetHLTResult("HLT_DoubleEle8_CaloIdL_TrkIdVL_HT160_v3") )        	return true;
	if ( GetHLTResult("HLT_DoubleEle8_CaloIdL_TrkIdVL_HT150_v1") )        	return true;
	if ( GetHLTResult("HLT_DoubleEle8_CaloIdL_TrkIdVL_HT150_v2") )        	return true;
	if ( GetHLTResult("HLT_DoubleEle8_CaloIdL_TrkIdVL_HT150_v3") )        	return true;
	if ( GetHLTResult("HLT_DOUBLEELE8_HT160_TIGHT") )        		return true;
	if ( GetHLTResult("HLT_DoubleEle8_CaloIdT_TrkIdVL_HT160_v1") )        	return true;
	if ( GetHLTResult("HLT_DoubleEle8_CaloIdT_TrkIdVL_HT160_v2") )        	return true;
	if ( GetHLTResult("HLT_DoubleEle8_CaloIdT_TrkIdVL_HT160_v3") )        	return true;
	if ( GetHLTResult("HLT_DoubleEle8_CaloIdT_TrkIdVL_HT150_v1") )        	return true;
	if ( GetHLTResult("HLT_DoubleEle8_CaloIdT_TrkIdVL_HT150_v2") )        	return true;
	if ( GetHLTResult("HLT_DoubleEle8_CaloIdT_TrkIdVL_HT150_v3") )        	return true;
  }
  if(iAnalysis==2) {//mt2
	return passAnyMT2Trigger();
  }
  return false;

}

//________________________________________________________________________________
const bool RunEfficiency::passMuTriggers(int iAnalysis) {
  if(iAnalysis==0) {//JZB
	if ( GetHLTResult("HLT_DoubleMu6_v1") )        return true;
	if ( GetHLTResult("HLT_DoubleMu6_v2") )        return true;
	if ( GetHLTResult("HLT_DoubleMu6_v3") )        return true;
	if ( GetHLTResult("HLT_DoubleMu7_v1") )        return true;
	if ( GetHLTResult("HLT_DoubleMu7_v2") )        return true;
	if ( GetHLTResult("HLT_DoubleMu7_v3") )        return true;
  }
  if(iAnalysis==1) {//same sign
	if ( GetHLTResult("HLT_MU8_JET40") )        return true;
	if ( GetHLTResult("HLT_Mu8_Jet40_v1") )        return true;
	if ( GetHLTResult("HLT_Mu8_Jet40_v2") )        return true;
	if ( GetHLTResult("HLT_Mu8_Jet40_v3") )        return true;
	if ( GetHLTResult("HLT_Mu8_Jet40_v4") )        return true;
	if ( GetHLTResult("HLT_Mu8_Jet40_v5") )        return true;
	if ( GetHLTResult("HLT_Mu8_Jet40_v6") )        return true;
	if ( GetHLTResult("HLT_DOUBLEMU7") )        return true;
	if ( GetHLTResult("HLT_DoubleMu6_v1") )        return true;
	if ( GetHLTResult("HLT_DoubleMu6_v2") )        return true;
	if ( GetHLTResult("HLT_DoubleMu6_v3") )        return true;
	if ( GetHLTResult("HLT_DoubleMu7_v1") )        return true;
	if ( GetHLTResult("HLT_DoubleMu7_v2") )        return true;
	if ( GetHLTResult("HLT_DoubleMu7_v3") )        return true;
	if ( GetHLTResult("HLT_MU13_MU8") )        return true;
	if ( GetHLTResult("HLT_Mu13_Mu8_v1") )        return true;
	if ( GetHLTResult("HLT_Mu13_Mu8_v2") )        return true;
	if ( GetHLTResult("HLT_MU17_ELE8") )        return true;
	if ( GetHLTResult("HLT_Mu17_Ele8_CaloIdL_v1") )        return true;
	if ( GetHLTResult("HLT_Mu17_Ele8_CaloIdL_v2") )        return true;
	if ( GetHLTResult("HLT_Mu17_Ele8_CaloIdL_v3") )        return true;
	if ( GetHLTResult("HLT_Mu17_Ele8_CaloIdL_v4") )        return true;
	if ( GetHLTResult("HLT_Mu17_Ele8_CaloIdL_v5") )        return true;
	if ( GetHLTResult("HLT_MU8_ELE17") )        return true;
	if ( GetHLTResult("HLT_Mu8_Ele17_CaloIdL_v1") )        return true;
	if ( GetHLTResult("HLT_Mu8_Ele17_CaloIdL_v2") )        return true;
	if ( GetHLTResult("HLT_Mu8_Ele17_CaloIdL_v3") )        return true;
	if ( GetHLTResult("HLT_Mu8_Ele17_CaloIdL_v4") )        return true;
	if ( GetHLTResult("HLT_Mu8_Ele17_CaloIdL_v5") )        return true;
	if ( GetHLTResult("HLT_DOUBLEMU3_HT160") )        return true;
	if ( GetHLTResult("HLT_DoubleMu3_HT160_v2") )        return true;
	if ( GetHLTResult("HLT_DoubleMu3_HT160_v3") )        return true;
	if ( GetHLTResult("HLT_DoubleMu3_HT150_v1") )        return true;
	if ( GetHLTResult("HLT_DoubleMu3_HT150_v2") )        return true;
	if ( GetHLTResult("HLT_DoubleMu3_HT150_v3") )        return true;
  }
  if(iAnalysis==2) {//mt2
	return passAnyMT2Trigger();
  }
  return false;
} 

//______________________________________________________________________________
const bool RunEfficiency::passEMuTriggers(int iAnalysis) {
  if(iAnalysis==0) {//JZB
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
  if(iAnalysis==1) {//same sign
	if ( GetHLTResult("HLT_MU3_ELE8_HT160") )        			return true;
	if ( GetHLTResult("HLT_Mu3_Ele8_CaloIdL_TrkIdVL_HT160_v1") )        	return true;
	if ( GetHLTResult("HLT_Mu3_Ele8_CaloIdL_TrkIdVL_HT160_v2") )        	return true;
	if ( GetHLTResult("HLT_Mu3_Ele8_CaloIdL_TrkIdVL_HT160_v3") )        	return true;
	if ( GetHLTResult("HLT_Mu3_Ele8_CaloIdL_TrkIdVL_HT150_v1") )        	return true;
	if ( GetHLTResult("HLT_Mu3_Ele8_CaloIdL_TrkIdVL_HT150_v2") )        	return true;
	if ( GetHLTResult("HLT_Mu3_Ele8_CaloIdL_TrkIdVL_HT150_v3") )        	return true;
	if ( GetHLTResult("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v1") )        return true;
	if ( GetHLTResult("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v2") )        return true;
	if ( GetHLTResult("HLT_MU3_ELE8_HT160_TIGHT") )        			return true;
	if ( GetHLTResult("HLT_Mu3_Ele8_CaloIdT_TrkIdVL_HT160_v1") )        	return true;
	if ( GetHLTResult("HLT_Mu3_Ele8_CaloIdT_TrkIdVL_HT160_v2") )        	return true;
	if ( GetHLTResult("HLT_Mu3_Ele8_CaloIdT_TrkIdVL_HT160_v3") )        	return true;
	if ( GetHLTResult("HLT_Mu3_Ele8_CaloIdT_TrkIdVL_HT150_v1") )        	return true;
	if ( GetHLTResult("HLT_Mu3_Ele8_CaloIdT_TrkIdVL_HT150_v2") )        	return true;
	if ( GetHLTResult("HLT_Mu3_Ele8_CaloIdT_TrkIdVL_HT150_v3") )        	return true;
  }
  if(iAnalysis==2) {//mt2
	return passAnyMT2Trigger();
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
	if(passElTriggers(0)) t_nEvent.passedee_triggersjzb = 1;
	if(passMuTriggers(0)) t_nEvent.passedmm_triggersjzb = 1;
	if(passEMuTriggers(0)) t_nEvent.passedem_triggersjzb = 1;

	if(passElTriggers(1)) t_nEvent.passedee_triggersss = 1;
	if(passMuTriggers(1)) t_nEvent.passedmm_triggersss = 1;
	if(passEMuTriggers(1)) t_nEvent.passedem_triggersss = 1;

	if(passElTriggers(2)) t_nEvent.passedee_triggersmt2 = 1;
	if(passMuTriggers(2)) t_nEvent.passedmm_triggersmt2 = 1;
	if(passEMuTriggers(2)) t_nEvent.passedem_triggersmt2 = 1;
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

      tmpLepton.tagjzb = 0;		tmpLepton.tagss = 0;		tmpLepton.tagmt2 = 0;
      tmpLepton.probeRecojzb = 0;	tmpLepton.probeRecoss = 0;	tmpLepton.probeRecomt2 = 0;
      tmpLepton.pprobeRecojzb = 0;	tmpLepton.pprobeRecoss = 0;	tmpLepton.pprobeRecomt2 = 0;
      tmpLepton.probeIsojzb = 0;	tmpLepton.probeIsoss = 0;	tmpLepton.probeIsomt2 = 0;
      tmpLepton.pprobeIsojzb = 0;	tmpLepton.pprobeIsoss = 0;	tmpLepton.pprobeIsomt2 = 0;
      tmpLepton.probeIDjzb = 0;		tmpLepton.probeIDss = 0;	tmpLepton.probeIDmt2 = 0;
      tmpLepton.pprobeIDjzb = 0;	tmpLepton.pprobeIDss = 0;	tmpLepton.pprobeIDmt2 = 0;

      if(MuPassingTag(0, muIndex)) tmpLepton.tagjzb = 1;
      if(MuPassingTag(1, muIndex)) tmpLepton.tagss = 1;
      if(MuPassingTag(2, muIndex)) tmpLepton.tagmt2 = 1;
      
      if(MuPassingRecoProbe(0, muIndex)) tmpLepton.probeRecojzb = 1;
      if(MuPassingRecoProbe(1, muIndex)) tmpLepton.probeRecoss = 1;
      if(MuPassingRecoProbe(2, muIndex)) tmpLepton.probeRecomt2 = 1;

      if(MuPassingRecoPProbe(0, muIndex)) tmpLepton.pprobeRecojzb = 1;
      if(MuPassingRecoPProbe(1, muIndex)) tmpLepton.pprobeRecoss = 1;
      if(MuPassingRecoPProbe(2, muIndex)) tmpLepton.pprobeRecomt2 = 1;

      if(MuPassingIDProbe(0, muIndex)) tmpLepton.probeIDjzb = 1;
      if(MuPassingIDProbe(1, muIndex)) tmpLepton.probeIDss = 1;
      if(MuPassingIDProbe(2, muIndex)) tmpLepton.probeIDmt2 = 1;

      if(MuPassingIDPProbe(0, muIndex)) tmpLepton.pprobeIDjzb = 1;
      if(MuPassingIDPProbe(1, muIndex)) tmpLepton.pprobeIDss = 1;
      if(MuPassingIDPProbe(2, muIndex)) tmpLepton.pprobeIDmt2 = 1;
      
      if(MuPassingIsoProbe(0, muIndex)) tmpLepton.probeIsojzb = 1;
      if(MuPassingIsoProbe(1, muIndex)) tmpLepton.probeIsoss = 1;
      if(MuPassingIsoProbe(2, muIndex)) tmpLepton.probeIsomt2 = 1;

      if(MuPassingIsoPProbe(0, muIndex, indexOfAssociatedPFMuon)) tmpLepton.pprobeIsojzb = 1;
      if(MuPassingIsoPProbe(1, muIndex, indexOfAssociatedPFMuon)) tmpLepton.pprobeIsoss = 1;
      if(MuPassingIsoPProbe(2, muIndex, indexOfAssociatedPFMuon)) tmpLepton.pprobeIsomt2 = 1;

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

      tmpLepton.tagjzb = 0;		tmpLepton.tagss = 0;		tmpLepton.tagmt2 = 0;
      tmpLepton.probeRecojzb = 0;	tmpLepton.probeRecoss = 0;	tmpLepton.probeRecomt2 = 0;
      tmpLepton.pprobeRecojzb = 0;	tmpLepton.pprobeRecoss = 0;	tmpLepton.pprobeRecomt2 = 0;
      tmpLepton.probeIsojzb = 0;	tmpLepton.probeIsoss = 0;	tmpLepton.probeIsomt2 = 0;
      tmpLepton.pprobeIsojzb = 0;	tmpLepton.pprobeIsoss = 0;	tmpLepton.pprobeIsomt2 = 0;
      tmpLepton.probeIDjzb = 0;		tmpLepton.probeIDss = 0;	tmpLepton.probeIDmt2 = 0;
      tmpLepton.pprobeIDjzb = 0;	tmpLepton.pprobeIDss = 0;	tmpLepton.pprobeIDmt2 = 0;

      if(ElPassingTag(0, elIndex)) tmpLepton.tagjzb = 1;
      if(ElPassingTag(1, elIndex)) tmpLepton.tagss = 1;
      if(ElPassingTag(2, elIndex)) tmpLepton.tagmt2 = 1;

      if(ElPassingRecoProbe(0, elIndex)) tmpLepton.probeRecojzb = 1;
      if(ElPassingRecoProbe(1, elIndex)) tmpLepton.probeRecoss = 1;
      if(ElPassingRecoProbe(2, elIndex)) tmpLepton.probeRecomt2 = 1;

      if(ElPassingRecoPProbe(0, elIndex)) tmpLepton.pprobeRecojzb = 1;
      if(ElPassingRecoPProbe(1, elIndex)) tmpLepton.pprobeRecoss = 1;
      if(ElPassingRecoPProbe(2, elIndex)) tmpLepton.pprobeRecomt2 = 1;

      if(ElPassingIDProbe(0, elIndex)) tmpLepton.probeIDjzb = 1;
      if(ElPassingIDProbe(1, elIndex)) tmpLepton.probeIDss = 1;
      if(ElPassingIDProbe(2, elIndex)) tmpLepton.probeIDmt2 = 1;

      if(ElPassingIDPProbe(0, elIndex)) tmpLepton.pprobeIDjzb = 1;
      if(ElPassingIDPProbe(1, elIndex)) tmpLepton.pprobeIDss = 1;
      if(ElPassingIDPProbe(2, elIndex)) tmpLepton.pprobeIDmt2 = 1;

      if(ElPassingIsoProbe(0, elIndex)) tmpLepton.probeIsojzb = 1;
      if(ElPassingIsoProbe(1, elIndex)) tmpLepton.probeIsoss = 1;
      if(ElPassingIsoProbe(2, elIndex)) tmpLepton.probeIsomt2 = 1;

      if(ElPassingIsoPProbe(0, elIndex, indexOfAssociatedPFEl)) tmpLepton.pprobeIsojzb = 1;
      if(ElPassingIsoPProbe(1, elIndex, indexOfAssociatedPFEl)) tmpLepton.pprobeIsoss = 1;
      if(ElPassingIsoPProbe(2, elIndex, indexOfAssociatedPFEl)) tmpLepton.pprobeIsomt2 = 1;

    leptons.push_back(tmpLepton);
  
  }
  

  vector<t_lepton> sortedGoodLeptons = sortLeptonsByPt(leptons);

  if(sortedGoodLeptons.size() < 2) {
    if( isMC ) t_myTree->Fill();
    return;
  }

//  if(sortedGoodLeptons.size() > 1) {
    
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
    if(!(sortedGoodLeptons[PosLepton1].p.Pt() > 5 && sortedGoodLeptons[PosLepton2].p.Pt() > 5 && 
       fabs(sortedGoodLeptons[PosLepton1].p.Eta())< 2.4 && fabs(sortedGoodLeptons[PosLepton2].p.Eta())< 2.4)) {
      	if( isMC ) t_myTree->Fill();
      	return;
    }

      t_nEvent.eta1 = sortedGoodLeptons[PosLepton1].p.Eta();
      t_nEvent.pt1 = sortedGoodLeptons[PosLepton1].p.Pt();
      t_nEvent.phi1 = sortedGoodLeptons[PosLepton1].p.Phi();
      t_nEvent.ch1 = sortedGoodLeptons[PosLepton1].charge;
      t_nEvent.id1 = sortedGoodLeptons[PosLepton1].type; //??????
        t_nEvent.tag1jzb = sortedGoodLeptons[PosLepton1].tagjzb; //??????
        t_nEvent.tag1ss = sortedGoodLeptons[PosLepton1].tagss; //??????
        t_nEvent.tag1mt2 = sortedGoodLeptons[PosLepton1].tagmt2; //??????
        t_nEvent.probeReco1jzb = sortedGoodLeptons[PosLepton1].probeRecojzb; //??????
        t_nEvent.probeReco1ss = sortedGoodLeptons[PosLepton1].probeRecoss; //??????
        t_nEvent.probeReco1mt2 = sortedGoodLeptons[PosLepton1].probeRecomt2; //??????
        t_nEvent.pprobeReco1jzb = sortedGoodLeptons[PosLepton1].pprobeRecojzb; //??????
        t_nEvent.pprobeReco1ss = sortedGoodLeptons[PosLepton1].pprobeRecoss; //??????
        t_nEvent.pprobeReco1mt2 = sortedGoodLeptons[PosLepton1].pprobeRecomt2; //??????
        t_nEvent.probeIso1jzb = sortedGoodLeptons[PosLepton1].probeIsojzb; //??????
        t_nEvent.probeIso1ss = sortedGoodLeptons[PosLepton1].probeIsoss; //??????
        t_nEvent.probeIso1mt2 = sortedGoodLeptons[PosLepton1].probeIsomt2; //??????
        t_nEvent.pprobeIso1jzb = sortedGoodLeptons[PosLepton1].pprobeIsojzb; //??????
        t_nEvent.pprobeIso1ss = sortedGoodLeptons[PosLepton1].pprobeIsoss; //??????
        t_nEvent.pprobeIso1mt2 = sortedGoodLeptons[PosLepton1].pprobeIsomt2; //??????
        t_nEvent.probeID1jzb = sortedGoodLeptons[PosLepton1].probeIDjzb; //??????
        t_nEvent.probeID1ss = sortedGoodLeptons[PosLepton1].probeIDss; //??????
        t_nEvent.probeID1mt2 = sortedGoodLeptons[PosLepton1].probeIDmt2; //??????
        t_nEvent.pprobeID1jzb = sortedGoodLeptons[PosLepton1].pprobeIDjzb; //??????
        t_nEvent.pprobeID1ss = sortedGoodLeptons[PosLepton1].pprobeIDss; //??????
        t_nEvent.pprobeID1mt2 = sortedGoodLeptons[PosLepton1].pprobeIDmt2; //??????
        t_nEvent.probe1jzb = t_nEvent.probeReco1jzb*t_nEvent.probeIso1jzb*t_nEvent.probeID1jzb;
        t_nEvent.probe1ss = t_nEvent.probeReco1ss*t_nEvent.probeIso1ss*t_nEvent.probeID1ss;
        t_nEvent.probe1mt2 = t_nEvent.probeReco1mt2*t_nEvent.probeIso1mt2*t_nEvent.probeID1mt2;
        t_nEvent.pprobe1jzb = t_nEvent.pprobeReco1jzb*t_nEvent.pprobeIso1jzb*t_nEvent.pprobeID1jzb;
        t_nEvent.pprobe1ss = t_nEvent.pprobeReco1ss*t_nEvent.pprobeIso1ss*t_nEvent.pprobeID1ss;
        t_nEvent.pprobe1mt2 = t_nEvent.pprobeReco1mt2*t_nEvent.pprobeIso1mt2*t_nEvent.pprobeID1mt2;

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
//      for(int n = 0; n < nAnalysis; ++n) {
        t_nEvent.tag2jzb = sortedGoodLeptons[PosLepton2].tagjzb; //??????
        t_nEvent.tag2ss = sortedGoodLeptons[PosLepton2].tagss; //??????
        t_nEvent.tag2mt2 = sortedGoodLeptons[PosLepton2].tagmt2; //??????
        t_nEvent.probeReco2jzb = sortedGoodLeptons[PosLepton2].probeRecojzb; //??????
        t_nEvent.probeReco2ss = sortedGoodLeptons[PosLepton2].probeRecoss; //??????
        t_nEvent.probeReco2mt2 = sortedGoodLeptons[PosLepton2].probeRecomt2; //??????
        t_nEvent.pprobeReco2jzb = sortedGoodLeptons[PosLepton2].pprobeRecojzb; //??????
        t_nEvent.pprobeReco2ss = sortedGoodLeptons[PosLepton2].pprobeRecoss; //??????
        t_nEvent.pprobeReco2mt2 = sortedGoodLeptons[PosLepton2].pprobeRecomt2; //??????
        t_nEvent.probeIso2jzb = sortedGoodLeptons[PosLepton2].probeIsojzb; //??????
        t_nEvent.probeIso2ss = sortedGoodLeptons[PosLepton2].probeIsoss; //??????
        t_nEvent.probeIso2mt2 = sortedGoodLeptons[PosLepton2].probeIsomt2; //??????
        t_nEvent.pprobeIso2jzb = sortedGoodLeptons[PosLepton2].pprobeIsojzb; //??????
        t_nEvent.pprobeIso2ss = sortedGoodLeptons[PosLepton2].pprobeIsoss; //??????
        t_nEvent.pprobeIso2mt2 = sortedGoodLeptons[PosLepton2].pprobeIsomt2; //??????
        t_nEvent.probeID2jzb = sortedGoodLeptons[PosLepton2].probeIDjzb; //??????
        t_nEvent.probeID2ss = sortedGoodLeptons[PosLepton2].probeIDss; //??????
        t_nEvent.probeID2mt2 = sortedGoodLeptons[PosLepton2].probeIDmt2; //??????
        t_nEvent.pprobeID2jzb = sortedGoodLeptons[PosLepton2].pprobeIDjzb; //??????
        t_nEvent.pprobeID2ss = sortedGoodLeptons[PosLepton2].pprobeIDss; //??????
        t_nEvent.pprobeID2mt2 = sortedGoodLeptons[PosLepton2].pprobeIDmt2; //??????
        t_nEvent.probe2jzb = t_nEvent.probeReco2jzb*t_nEvent.probeIso2jzb*t_nEvent.probeID2jzb;
        t_nEvent.probe2ss = t_nEvent.probeReco2ss*t_nEvent.probeIso2ss*t_nEvent.probeID2ss;
        t_nEvent.probe2mt2 = t_nEvent.probeReco2mt2*t_nEvent.probeIso2mt2*t_nEvent.probeID2mt2;
        t_nEvent.pprobe2jzb = t_nEvent.pprobeReco2jzb*t_nEvent.pprobeIso2jzb*t_nEvent.pprobeID2jzb;
        t_nEvent.pprobe2ss = t_nEvent.pprobeReco2ss*t_nEvent.pprobeIso2ss*t_nEvent.pprobeID2ss;
        t_nEvent.pprobe2mt2 = t_nEvent.pprobeReco2mt2*t_nEvent.pprobeIso2mt2*t_nEvent.pprobeID2mt2;
//      }    
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
/*
    } else {
      
      if( isMC ) t_myTree->Fill();
      return;
      
    }*/

    
    float jetThreshold[nAnalysis] = {30, 40, 20};
//    for(int n = 0; n < nAnalysis; n++) { 
      t_nEvent.pfJetGoodNumjzb=0;
      t_nEvent.pfJetGoodNumss=0;
      t_nEvent.pfJetGoodNummt2=0;
//      vector<t_lepton> pfGoodJets;
      float closest1jzb=1000, closest2jzb=1000, closestnjzb=1000; //the dr distance between the leptons and jets
      float closest1ss=1000,  closest2ss=1000,  closestnss=1000; //the dr distance between the leptons and jets
      float closest1mt2=1000, closest2mt2=1000, closestnmt2=1000; //the dr distance between the leptons and jets
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
	bool isJetIDjzb = IsGoodBasicPFJet(i, false); //JZB
	bool isJetIDss  = IsGoodBasicPFJet(i, false); //SS
	bool isJetIDmt2 = IsGoodBasicPFJetPAT3(i, 20.0, 3.0); //MT2
	
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
	
	t_nEvent.pfHTjzb    += jpt;
	t_nEvent.pfHTss    += jpt;
	t_nEvent.pfHTmt2    += jpt;
	
	// Keep central jets
	if ( !(fabs(jeta)<2.6 ) ) continue;
	
	
	t_nEvent.pfGoodHTjzb += jpt;
	t_nEvent.pfGoodHTss += jpt;
	t_nEvent.pfGoodHTmt2 += jpt;

	sumOfPFJets += aJet;
	
/*	t_lepton tmpLepton;
	tmpLepton.p = aJet;
	tmpLepton.charge = 0;
	tmpLepton.index = i;
	tmpLepton.type = -1;*/
//	pfGoodJets.push_back(tmpLepton);
	
	if ( jpt > jetThreshold[0] ) {//JZB
	  t_nEvent.pfTightHTjzb += jpt;
	  if(isJetIDjzb>0) t_nEvent.pfJetGoodNumIDjzb++;
	  t_nEvent.pfJetGoodNumjzb++;
	  float dr1 = aJet.DeltaR(sortedGoodLeptons[PosLepton1].p);
	  float dr2 = aJet.DeltaR(sortedGoodLeptons[PosLepton2].p);
	  float drn = aJet.DeltaR(sortedGoodLeptons[negativePosition].p);
	  if(dr1<closest1jzb) closest1jzb = dr1;
	  if(dr2<closest2jzb) closest2jzb = dr2;
	  if(drn<closestnjzb) closestnjzb = drn;
	}
	if ( jpt > jetThreshold[1] ) {//SS
	  t_nEvent.pfTightHTss += jpt;
	  if(isJetIDss>0) t_nEvent.pfJetGoodNumIDss++;
	  t_nEvent.pfJetGoodNumss++;
	  float dr1 = aJet.DeltaR(sortedGoodLeptons[PosLepton1].p);
	  float dr2 = aJet.DeltaR(sortedGoodLeptons[PosLepton2].p);
	  float drn = aJet.DeltaR(sortedGoodLeptons[negativePosition].p);
	  if(dr1<closest1ss) closest1ss = dr1;
	  if(dr2<closest2ss) closest2ss = dr2;
	  if(drn<closestnss) closestnss = drn;
	}
	if ( jpt > jetThreshold[2] ) {//MT2
	  t_nEvent.pfTightHTmt2 += jpt;
	  if(isJetIDmt2>0) t_nEvent.pfJetGoodNumIDmt2++;
	  t_nEvent.pfJetGoodNummt2++;
	  float dr1 = aJet.DeltaR(sortedGoodLeptons[PosLepton1].p);
	  float dr2 = aJet.DeltaR(sortedGoodLeptons[PosLepton2].p);
	  float drn = aJet.DeltaR(sortedGoodLeptons[negativePosition].p);
	  if(dr1<closest1mt2) closest1mt2 = dr1;
	  if(dr2<closest2mt2) closest2mt2 = dr2;
	  if(drn<closestnmt2) closestnmt2 = drn;
	}
	
      }
      
      t_nEvent.pfMETjzb = fTR->PFMET;
      t_nEvent.pfMETss = fTR->PFMET;
      t_nEvent.pfMETmt2 = fTR->PFMET;

      t_nEvent.dr1jzb = closest1jzb;
      t_nEvent.dr1ss = closest1ss;
      t_nEvent.dr1mt2 = closest1mt2;
      t_nEvent.dr2jzb = closest2jzb;
      t_nEvent.dr2ss = closest2ss;
      t_nEvent.dr2mt2 = closest2mt2;
      t_nEvent.drnjzb = closestnjzb;
      t_nEvent.drnss = closestnss;
      t_nEvent.drnmt2 = closestnmt2;
//    }




      t_myTree->Fill();	
//  } else {
  
//    if( isMC ) t_myTree->Fill();
   

//  }

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


//_____________________________________________________________________________
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

  float closest1 = 1000, closest2 = 1000, dr1, dr2;
  if(sortedGLeptons.size()>1) { 
    // Number of good jets
    t_nEvent.genNjets = 0;
    for ( int jIndex=0; jIndex<fTR->NGenJets; ++jIndex) {
      if ( fTR->GenJetPt[jIndex]<minJPt ) continue;
      if ( fabs(fTR->GenJetEta[jIndex])>maxJEta ) continue;
      ++t_nEvent.genNjets;
      TLorentzVector tmpVector;
      tmpVector.SetPtEtaPhiE(fTR->GenJetPt[jIndex],fTR->GenJetEta[jIndex], fTR->GenJetPhi[jIndex], fTR->GenJetE[jIndex]);
      dr1 = tmpVector.DeltaR(sortedGLeptons[i1].p);
      dr2 = tmpVector.DeltaR(sortedGLeptons[i2].p);
      if(dr1<closest1) closest1 = dr1;
      if(dr2<closest2) closest2 = dr2;
    }
    dr1 = closest1;
    dr2 = closest2;
  }
  

  
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
