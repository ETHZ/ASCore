#include "helper/Utilities.hh"
#include "MassAnalysis.hh"
#include "TLorentzVector.h"
#include <sstream>

using namespace std;

MassAnalysis::MassAnalysis(TreeReader *tr) : MultiplicityAnalysisBase(tr){
	Util::SetStyle();	
}

MassAnalysis::~MassAnalysis(){
}

void MassAnalysis::Begin(){
	// Define the output file of histograms
	const char* filename = "Mass_histos.root";
	fHistFile = new TFile(fOutputDir + TString(filename), "RECREATE");

	fETHnt = new ETHnt();
	// book tree
	BookTree();

		
	fMT2_histos_step = 50;
	fMT2_histos_number = 6;
	
	// Define the histograms
	fHMT2_vs_M_OSDiLept = new TH2D("MT2_vs_M_OSDiLept"   , "MT2_vs_M_OSDiLept"       , 200, 0, 200, 100, 0, 500);
	fHAlpahT_DiJet      = new TH1D("AlpahT_DiJet"        , "alpha_T for dijets"      , 300, 0., 2.);
	fHAlpahT_PseudoJet  = new TH1D("AlpahT_PseudoJet"    , "alpha_T for pseudojets"  , 300, 0., 2.);
	
	fHMCT_DiJet                  = new TH1D("MCT_DiJet"                , "di-Jets MCT"                ,  60, 0., 500.);
	fHMCT_DiBJet                 = new TH1D("MCT_DiBJet"               , "di-B-Jets MCT"              ,  50, 0., 500.);
	fHMCT_PseudoJet              = new TH1D("MCT_PseudoJet"            , "pseudo jet MCT"             , 150, 0., 1500.);
	fHMCT_PseudoJetWithLeptons   = new TH1D("MCT_PseudoJetWithLeptons" , "pseudo jet MCT with leptons",  50, 0., 1000.);
	fHMCT_PseudoJetNoLeptons     = new TH1D("MCT_PseudoJetNoLeptons"   , "pseudo jet MCT lepton veto ", 150, 0., 1500.);
	fHMCT_PseudoJetWithB         = new TH1D("MCT_PseudoJetWithB"       , "pseudo jet MCT with b"      , 100, 0., 1500.);
	fHMCT_OSee                   = new TH1D("MCT_OSee"           , "OS ee MCT"                ,  40, 0., 400.);
	fHMCT_OSmumu                 = new TH1D("MCT_OSmumu"         , "OS mumu MCT"              ,  40, 0., 400.);
	fHMCT_OSemu                  = new TH1D("MCT_OSemu"          , "OS emu MCT"               ,  40, 0., 400.);
	                                                                                 
	fHMCTcorr_DiJet                  = new TH1D("MCTcorr_DiJet"                , "di-Jets MCT"                ,  60, 0., 400.);
	fHMCTcorr_DiBJet                 = new TH1D("MCTcorr_DiBJet"               , "di-B-Jets MCT"              ,  30, 0., 300.);
	fHMCTcorr_PseudoJet              = new TH1D("MCTcorr_PseudoJet"            , "pseudo jet MCT"             , 200, 0., 1500.);
	fHMCTcorr_PseudoJetWithLeptons   = new TH1D("MCTcorr_PseudoJetWithLeptons" , "pseudo jet MCT with leptons",  50, 0., 1000.);
	fHMCTcorr_PseudoJetNoLeptons     = new TH1D("MCTcorr_PseudoJetNoLeptons"   , "pseudo jet MCT lepton veto ", 200, 0., 1500.);
	fHMCTcorr_PseudoJetWithB         = new TH1D("MCTcorr_PseudoJetWithB"       , "pseudo jet MCT with b"      , 200, 0., 1500.);
	fHMCTcorr_OSee                   = new TH1D("MCTcorr_OSee"           , "OS ee MCT"                ,  20, 0., 200.);
	fHMCTcorr_OSmumu                 = new TH1D("MCTcorr_OSmumu"         , "OS mumu MCT"              ,  20, 0., 200.);
	fHMCTcorr_OSemu                  = new TH1D("MCTcorr_OSemu"          , "OS emu MCT"               ,  20, 0., 200.);


	fHMCTperp_OSee      = new TH1D("MCTperp_OSee"        , "OS ee MCT perp"          ,  20, 0., 200.);	
	fHMCTperp_OSmumu    = new TH1D("MCTperp_OSmumu"      , "OS mumu MCT perp"        ,  20, 0., 200.);	
	fHMCTperp_OSemu     = new TH1D("MCTperp_OSemu"       , "OS emu MCT perp"         ,  20, 0., 200.);	    
	                                             
	fHMCT_TTbar         = new TH1D("MCT_TTbar"           , "MCT for TTbar on GenLevel"      , 50, 0., 200.);
	fHMCTperp_TTbar     = new TH1D("MCTperp_TTbar"       , "TTbar MCT perp"                 , 50, 0., 200.);
	fHMCTcorr_TTbar     = new TH1D("MCTcorr_TTbar"       , "MCTcorr for TTbar on GenLevel"  , 50, 0., 200.);
	fHMT2_TTbar         = new TH1D("MT2_TTbar"           , "MT2 for TTbar with Gen Info"    , 50, 0., 200.);
	fHMT2perp_TTbar     = new TH1D("MT2perp_TTbar"       , "MT2perp for TTbar with Gen Info", 50, 0., 200.);
	
	fHInvMassOSee      = new TH1D("InvMass_OSee"  ,  " Inv Mass of OS di-leptons ee"   , 30, 0, 300);
	fHInvMassOSmumu    = new TH1D("InvMass_OSmumu",  " Inv Mass of OS di-leptons mumu" , 30, 0, 300);
	fHInvMassOSemu     = new TH1D("InvMass_OSemu" ,  " Inv Mass of OS di-leptons emu"  , 30, 0, 300);
	fHInvMassOSll      = new TH1D("InvMass_OSll"  ,  " Inv Mass of OS di-leptons ll"   , 30, 0, 300);
	fHInvMassSSee      = new TH1D("InvMass_SSee"  ,  " Inv Mass of SS di-leptons ee"   , 30, 0, 300);
	fHInvMassSSmumu    = new TH1D("InvMass_SSmumu",  " Inv Mass of SS di-leptons mumu" , 30, 0, 300);
	fHInvMassSSemu     = new TH1D("InvMass_SSemu" ,  " Inv Mass of SS di-leptons emu"  , 30, 0, 300);
	fHInvMassSSll      = new TH1D("InvMass_SSll"  ,  " Inv Mass of SS di-leptons ll"   , 30, 0, 300);
	fHInvMassDijet     = new TH1D("InvMassDijet"  ,  " Inv Mass of cleaned Di-Jets"    ,100, 0, 2500);
	fHInvMassDiBjet    = new TH1D("InvMassDiBjet" ,  " Inv Mass of cleaned Di-BJets"   , 50, 0, 1000);
	fHInvMassdiBHemi   = new TH1D("InvMassdiBHemi",  " Inv Mass of two b jets in same hemosphere", 20, 0, 300);

	fHDPhiJ1vsDPhiJ2         = new TH2D("DPhiJ1vsDPhiJ2"      ,  "DPhiJ1 vs DPhiJ2"         , 100., 0., 3.2, 100., 0., 3.2);
	fHDPhiJ1vsDPhiJ2_clean   = new TH2D("DPhiJ1vsDPhiJ2_clean",  "DPhiJ1 vs DPhiJ2"         , 100., 0., 3.2, 100., 0., 3.2);
	fHMHT                    = new TH1D("MHT"      , "MHT of jets only"                     , 100., 0., 700.);
	fHMHT_clean              = new TH1D("MHT_clean", "MHT of jets only"                     , 100., 0., 700.);
	fHPFMET_clean            = new TH1D("PFMET_clean"         ,  "PFMET"                    ,      100, 0., 500.);
	fHPFMET                  = new TH1D("PFMET"               ,  "PFMET"                    ,      100, 0., 500.);
	fHJpt_clean              = new TH1D("JPt_clean"           ,  "Jet Pt, all selected jets",      200, 0., 1000);  
	fHJpt                    = new TH1D("JPt"                 ,  "Jet Pt, all selected jets",      200, 0., 1000);  
	fHJEta_clean             = new TH1D("JEta_clean"          ,  "Jet Eta, all selected jets",     100, -3,   3);
	fHJEta                   = new TH1D("JEta"                ,  "Jet Eta, all selected jets",     100, -3,   3);
	fHBJpt                   = new TH1D("JBPt"                ,  "BJet Pt, all selected jets",     200, 0., 1000);
	fHBJEta                  = new TH1D("JBEta"               ,  "BJet Eta, all selected jets",    100, -3,   3);
	fHNJets                  = new TH1D("NJets",                 "NJets",                          15, 0., 15. );
	fHNJets_clean            = new TH1D("NJets_clean",           "NJets",                          15, 0., 15. );
	fHLeptConfig             = new TH1D("LeptConfig"      ,      "LeptConfig",                     12, 0., 12. );
	fHLeptConfig_clean       = new TH1D("LeptConfig_clean",      "LeptConfig",                     12, 0., 12. );
	fHPseudoJetMT2AxisdPhi   = new TH1D("PseudoJetMT2AxisdPhi",  "PseudoJetMT2AxisdPhi",           100, 0.,3.2 );
	fHMPT                    = new TH1D("MPT",                   "MPT",                            100.,0.,500 );
	fHMPT_selected           = new TH1D("MPTselected",           "MPT",                            100.,0.,500 );

	fHVectorSumPt            = new TH1D("VectorSumPt"         ,  "norm of vector sum of Pt of selected objets + MET", 100, 0., 500.);
	fHVectorSumPt_clean      = new TH1D("VectorSumPt_clean"   ,  "norm of vector sum of Pt of selected objets + MET", 100, 0., 500.);
	fHVectorSumPtvsDiJetMT2     = new TH2D("VectorSumPtvsDiJetMT2", " VectorSumPT of Dijet MT2 M=0 ", 100, 0., 400., 100, 0., 300);
	fHPseudoJetMT2vsVectorSumPt = new TH2D("PseudoJetMT2vsVectorSumPt", " VectorSumPT of Pseudojet MT2 M=0 ", 100, 0., 300., 100, 0., 400);
	fHPseudoJetMT2vsMETsign     = new TH2D("PseudoJetMT2vsMETsign", " METsign of Pseudojet MT2 M=0 ", 100, 0., 300., 60, 0., 20);
	fHPseudoJetMT2vsMET         = new TH2D("PseudoJetMT2vsMET", "MT2 pseudojet vs MET", 100, 0., 300., 100, 0., 400);
	fHPseudoJetMT2vsLeadingJEta = new TH2D("PseudoJetMT2vsLeadingJEta", "MT2 pseudojet vs leading JEta", 100, 0., 300., 100, -3., 3.);
	fHPseudoJetMT2vsAlphaT      = new TH2D("PseudoJetMT2vsAlphaT", "MT2 pseudojet vs AlphaT", 100, 0., 300., 100, 0., 2.);
	fHPseudoJetMT2vsMHT         = new TH2D("PseudoJetMT2vsMHT", "MT2 pseudojet vs jet MHT", 100, 0., 300., 100, 0., 700.);	

	fHMT_single_e         = new TH1D("MT_single_e "  ,  " transverse mass of el, MET "    ,100, 0, 200);
	fHMT_single_mu	      = new TH1D("MT_single_mu"  ,  " transverse mass of mu, MET "    ,100, 0, 200);
	fHMT_single_e_nojets  = new TH1D("MT_single_e_nojets "  ,  " transverse mass of el, MET, no jets in event "    ,100, 0, 200);
	fHMT_single_mu_nojets = new TH1D("MT_single_mu_nojets"  ,  " transverse mass of mu, MET, no jets in event "    ,100, 0, 200);

	fHPFandCalo_deltaR   = new TH1D("PFandCalo_deltaR", "", 200, 0., 10);
	fHPFJ1J2DeltaR       = new TH1D("PFJ1J2DeltaR", "", 200, 0, 10. );

	for(int i=0; i<fMT2_histos_number; ++i){
		std::stringstream out;
		int mass = i*fMT2_histos_step;
		out << mass ;
		TString Mass = (TString) out.str();
				
		fHMT2_diBjet[i]                = new TH1D("MT2_diBjet_M"+Mass,               " ", 50, 0., 300.);
		fHMT2_dijet[i]                 = new TH1D("MT2_dijet_M"+Mass,                " ", 80, 0., 500.);
		fHMT2_diBjetdiLept[i]          = new TH1D("MT2_diBjetdiLept_M"+Mass,         " ",100, 0., 500.);
		fHMT2_pseudojet[i]             = new TH1D("MT2_pseudojet_M"+Mass,            " ",100, 0., 500.);
		fHMT2_PseudoJetWithB[i]        = new TH1D("MT2_PseudoJetsWithB_M"+Mass,      " ",100, 0., 500.);     
		fHMT2_PseudoJetWithLeptons[i]  = new TH1D("MT2_PseudoJetsWithLeptons_M"+Mass," ", 50, 0., 500.);
		fHMT2_PseudoJetNoLeptons[i]    = new TH1D("MT2_PseudoJetsNoLeptons_M"+Mass,  " ",100, 0., 500.);


		fHMT2_OSll[i]                 = new TH1D("MT2_OSll_M"+Mass,               " ", 30, 0., 300.);
		fHMT2_OSee[i]                 = new TH1D("MT2_OSee_M"+Mass,               " ", 30, 0., 300.);
		fHMT2_OSmumu[i]               = new TH1D("MT2_OSmumu_M"+Mass,             " ", 30, 0., 300.);
		fHMT2_OSemu[i]                = new TH1D("MT2_OSemu_M"+Mass,              " ", 30, 0., 300.);
		fHMT2_OSllminusemu[i]         = new TH1D("MT2_OSllminusemu_M"+Mass,       " ", 30, 0., 300.);
		
		fHMT2_SSll[i]                 = new TH1D("MT2_SSll_M"+Mass,               " ", 10, 0., 200.);
		
		fHMT2perp_OSee[i]             = new TH1D("MT2perp_OSee_M"+Mass,           " ", 20, 0., 200.);
		fHMT2perp_OSmumu[i]           = new TH1D("MT2perp_OSmumu_M"+Mass,         " ", 20, 0., 200.);
		fHMT2perp_OSemu[i]            = new TH1D("MT2perp_OSemu_M"+Mass,          " ", 20, 0., 200.);
	}
}

void MassAnalysis::Analyze(){	

	// ---------------------------------------------------
	// Initialize fElecs, fJets, fBJets, fMuons, fLeptConfig 
	InitializeEvent();
	// ----------------------------------------------------
	
	// --------------------------------------------------------------------
	// check if event passes selection cuts and trigger requirements
	if(! IsSelectedEvent()){return;}
	// --------------------------------------------------------------------
		
	// reset tree variables
	ResetTree();
	
	// Controll Plots
	ControlPlots();

	// Single Leptons (e, mu)
	SingleLeptonMasses();

	// Mass distributions for OS dileptons
	DiLeptonMasses();

	// Masses for di-jets
	DiJetMasses();
	// Masses for di-bjets

	DiBJetMasses();
	
	// Masses for at least three 3 jets
	MultiBMasses();

	// Masses for pseudojets
	PseudoJetMasses();

	// Masses for TTbar on GenLevel
	MassesForTTbar();

	// Fill Tree !! has to be called at the very end !!
	MassAnalysis::FillTree();

}

// ***********************************************************************************************
void MassAnalysis::BookTree(){
	           
	fETHtree = new TTree("ETHtree", "ETHtree");
	fETHtree->Branch("ETHnt" , "ETHnt" , &fETHnt);

	fMassTree = new TTree("MassAnalysis", "MassAnalysisTree");
	
	fMassTree->Branch("Run"                   ,&fTrun,                    "Run/I");         
	fMassTree->Branch("Event"                 ,&fTevent,                  "Event/I");       
	fMassTree->Branch("LumiSection"           ,&fTlumisec,                "LumiSection/I");
	fMassTree->Branch("HCALNoiseFlag"         ,&fThcalnoiseflag,          "HCALNoiseFlag/I");
	fMassTree->Branch("Weight"                ,&fWeight,                  "Weight/F");
	fMassTree->Branch("NJets"                 ,&fTnjets,                  "NJets/I");       
	fMassTree->Branch("JPt"                   ,&fTjpt,                    "JPt[NJets]/D");       
	fMassTree->Branch("DPhiJetsMet"           ,&fTjetmetdphi,             "DPhiJetsMet[NJets]/D");       
	fMassTree->Branch("NElecs"                ,&fTnelecs,                 "NElecs/I");       
	fMassTree->Branch("NMuons"                ,&fTnmuons,                 "NMuons/I");       
	fMassTree->Branch("LeptConfig"            ,&fTleptconfig,             "LeptConfig/I");       
	fMassTree->Branch("IsCleanMultiJetEvent"  ,&fTisCleanMultiJetEvent,   "IsCleanMultiJetEvent/I");       
	fMassTree->Branch("IsCleanJetEvent"       ,&fTisCleanJetEvent,        "IsCleanJetEvent/I");       
	fMassTree->Branch("R12R21"                ,&fTr12r21,                 "R12R21/I");
	fMassTree->Branch("NJetsPt50Eta25"        ,&fTnJetsPt50Eta25,         "NJetsPt50Eta25/I");
	fMassTree->Branch("PseudoJetMT2"          ,&fTpseudoJetMT2,           "PseudoJetMT2/D");
	fMassTree->Branch("PseudoJetMT2massive"   ,&fTpseudoJetMT2massive,    "PseudoJetMT2massive/D");
	fMassTree->Branch("PseudoJetMT2simple"    ,&fTpseudoJetMT2simple,     "PseudoJetMT2simple/D");
	fMassTree->Branch("PseudoJetMCT"          ,&fTpseudoJetMCT,           "PseudoJetMCT/D");
	fMassTree->Branch("PseudoJetMCTMinusMass" ,&fTpseudoJetMCTMinusMass,  "PseudoJetMCTMinusMass/D");
	fMassTree->Branch("PseudoJetMCTMassless"  ,&fTpseudoJetMCTmassless ,  "PseudoJetMCTMassless/D");
	fMassTree->Branch("PseudoJet1Pt"          ,&fTpseudoJet1Pt,           "PseudoJet1Pt/D");
	fMassTree->Branch("PseudoJet2Pt"          ,&fTpseudoJet2Pt,           "PseudoJet2Pt/D");
	fMassTree->Branch("PseudojetAlphaT"       ,&fTpseudojetAlphaT,        "PseudojetAlphaT/D");
	fMassTree->Branch("Vectorsumpt"           ,&fTvectorsumpt,            "Vectorsumpt/D");
	fMassTree->Branch("PFMET"                 ,&fTpfmet,                  "PFMET/D");
	fMassTree->Branch("PFMETsign"             ,&fTpfmetsign,              "PFMETsign/D");
	fMassTree->Branch("LeadingJetEta"         ,&fTleadingJetEta,          "LeadingJetEta/D");
	fMassTree->Branch("DPhiJ1Met"             ,&fTdPhiJ1MET,              "DPhiJ1Met/D");
	fMassTree->Branch("DPhiJ2Met"             ,&fTdPhiJ2MET,              "DPhiJ2Met/D");
	fMassTree->Branch("PseudoJetMT2AxisdPhi"  ,&fTPseudoJetMT2AxisdPhi,   "PseudoJetMT2AxisdPhi/D");
	fMassTree->Branch("R1221min"              ,&fTr1221min,               "R1221min/D");
	fMassTree->Branch("MPT_sel"               ,&fTmpt_sel,                "MPT_sel/D");
	fMassTree->Branch("MPT"                   ,&fTmpt,                    "MPT/D");
	fMassTree->Branch("MHT"                   ,&fTmHT,                    "MHT/D");
	fMassTree->Branch("HT"                    ,&fThT,                     "HT/D");
	fMassTree->Branch("DPhiMhtMpt"            ,&fTdPhiMhtMpt,             "DPhiMhtMpt/D");
}

void MassAnalysis::ResetTree(){
        fETHnt->Reset();

	fTrun                  =-1;
	fTevent                =-1;
	fTlumisec              =-1;
	fThcalnoiseflag        =-1;
	fTweight               =-99999.99;
	fTnjets                =-1;
	fTnelecs               =-1;
	fTnmuons               =-1;
	fTleptconfig           =-1;
	fTisCleanMultiJetEvent =-1;
	fTisCleanJetEvent      =-1;
	fTr12r21               =-1;
	fTnJetsPt50Eta25       =-1;
	fTpseudoJetMT2         =-99999.99;
	fTpseudoJetMT2massive  =-99999.99;
	fTpseudoJetMT2simple   =-99999.99;
	fTpseudoJetMCT         =-99999.99;
	fTpseudoJetMCTMinusMass=-99999.99;
	fTpseudoJetMCTmassless =-99999.99;
	fTpseudoJet1Pt         =-99999.99;
	fTpseudoJet2Pt         =-99999.99;
	fTpseudojetAlphaT      =-99999.99;
	fTvectorsumpt          =-99999.99;
	fTmHT                  =-99999.99;
	fTmpt_sel              =-99999.99;
        fTmpt                  =-99999.99;	
	fThT                   =-99999.99;
	fTpfmet                =-99999.99;
	fTpfmetsign            =-99999.99;
	fTleadingJetEta        =-99999.99;
	fTdPhiJ1MET            =-99999.99;
	fTdPhiJ2MET            =-99999.99;
	fTPseudoJetMT2AxisdPhi =-99999.99;
	fTr1221min             =-99999.99;
	fTdPhiMhtMpt           =-99999.99;
	
	for(int i=0; i<gMaxnjets; ++i){
		fTjpt[i]       =-99999.99;
		fTjetmetdphi[i]=-99999.99;
	}

}

void MassAnalysis::FillTree(){
	
	fTrun                  = fTR->Run;
	fTevent                = fTR->Event;
	fTlumisec              = fTR->LumiSection;
	fThcalnoiseflag        = fTR->HBHENoiseFlag;
	fTweight               = fWeight;


	fTnjets                = fJets.size();
	fTnelecs               = fElecs.size();
	fTnmuons               = fMuons.size();
	fTleptconfig           = (int) fLeptConfig;
	fTisCleanMultiJetEvent = (int) fIsCleanMultiJetEvent;
	fTisCleanJetEvent      = (int) fIsCleanJetEvent;
	fTr12r21               = (int) fR12R21;
	fTnJetsPt50Eta25       = fNJetsPt50Eta25;
	fTvectorsumpt          = fVectorSumPt;
	fTmHT                  = fMHT;
	fThT                   = fHT;
	fTpfmet                = fTR->PFMET;
	fTpfmetsign            = (fTR->PFMET)/sqrt(fTR->SumEt);

	if(fJets.size()>1 && fCut_PFMET_min > 0) fTr1221min  = (fR12 <= fR21) ? fR12 : fR21;
	if(fJets.size()>0) fTleadingJetEta   = fTR-> PFJEta[fJets[0]];
	if(fJets.size()>0 && fCut_PFMET_min > 0) fTdPhiJ1MET = fDeltaPhi1; 
	if(fJets.size()>0 && fCut_PFMET_min > 0) fTdPhiJ2MET = fDeltaPhi2; 

	// Tracks
	TVector3 tracks(0.,0.,0.);
	TVector3 track(0.,0.,0.);
	for(int i=0; i< fTR->NTracks; ++i){
		track.SetPtEtaPhi(fabs(fTR->TrkPt[i]), fTR->TrkEta[i], fTR->TrkPhi[i]);
		tracks += track;
	}	
	fTmpt_sel= tracks.Pt();
	fTmpt    = sqrt((fTR->TrkPtSumx)*(fTR->TrkPtSumx) + (fTR->TrkPtSumy)*(fTR->TrkPtSumy));
	
	double MHTphi=fMHTphi;
	double MPTphi=tracks.Phi();
	fTdPhiMhtMpt = fabs(Util::DeltaPhi(MHTphi, MPTphi));

	// control variables
	// jpt
	for(int i=0; i<fJets.size(); ++i){
		fTjpt[i]       =fTR->PFJPt[fJets[i]];
		fTjetmetdphi[i]=Util::DeltaPhi(fTR->PFJPhi[fJets[i]], fTR->PFMETphi); 
	}

	fETHnt->misc.Run                   = fTR->Run;
	fETHnt->misc.Event		   = fTR->Event;
	fETHnt->misc.LumiSection	   = fTR->LumiSection;
	fETHnt->misc.Weight		   = fWeight;
	fETHnt->misc.LeptConfig		   = (int) fLeptConfig;
	fETHnt->misc.IsCleanMultiJetEvent  = fIsCleanMultiJetEvent;
	fETHnt->misc.IsCleanJetEvent	   = fIsCleanJetEvent;
	fETHnt->misc.R12R21		   = (int) fR12R21;
	fETHnt->misc.NJetsPt50Eta25	   = fNJetsPt50Eta25;
	fETHnt->misc.Vectorsumpt	   = fVectorSumPt;
	fETHnt->misc.MHT                   = fMHT;
	fETHnt->misc.HT	                   = fHT;
	fETHnt->misc.PFMET	  	   = fTR->PFMET;
	fETHnt->misc.PFMETphi	  	   = fTR->PFMETphi;
	fETHnt->misc.PFMETsign	           = (fTR->PFMET)/sqrt(fTR->SumEt);
	fETHnt->misc.MPT_sel	           = tracks.Pt();
	fETHnt->misc.MPT	           = sqrt((fTR->TrkPtSumx)*(fTR->TrkPtSumx) + (fTR->TrkPtSumy)*(fTR->TrkPtSumy));
	fETHnt->misc.DPhiMhtMpt            = fabs(Util::DeltaPhi(MHTphi, MPTphi));
	if(fJets.size()>1 && fCut_PFMET_min > 0) fETHnt->misc.R1221min        = (fR12 <= fR21) ? fR12 : fR21;
	if(fJets.size()>0)                       fETHnt->misc.LeadingJetEta   = fTR-> PFJEta[fJets[0]];
	if(fJets.size()>0 && fCut_PFMET_min > 0) fETHnt->misc.DPhiJ1Met       = fDeltaPhi1; 
	if(fJets.size()>0 && fCut_PFMET_min > 0) fETHnt->misc.DPhiJ2Met       = fDeltaPhi2; 


	fETHnt->SetNJets ((Int_t)fJets .size());
	fETHnt->SetNEles ((Int_t)fElecs.size());
	fETHnt->SetNMuons((Int_t)fMuons.size());

	// Fill jets 4-momenta
	for(int i=0; i<fJets.size(); ++i) {
	  fETHnt->jet[i].lv.SetPxPyPzE( fTR->PFJPx[fJets[i]], fTR->PFJPy[fJets[i]], 
					fTR->PFJPz[fJets[i]], fTR->PFJE [fJets[i]]); //SetLV(GetJet4Momenta(fJets[i]));
	  // b-tag info now should be available
	  // fETHnt->jet[i].bTagProbTCHE  =  fTR->PFJbTagProbTkCntHighEff [fJets[i]];
	  // fETHnt->jet[i].bTagProbTCHE  =  fTR->PFJbTagProbTkCntHighPur [fJets[i]];
	  // fETHnt->jet[i].bTagProbSSVHE =  fTR->PFJbTagProbSimpSVHighEff[fJets[i]];
	  // fETHnt->jet[i].bTagProbSSVHE =  fTR->PFJbTagProbSimpSVHighPur[fJets[i]];
	}
	// Fill leptons 4-momenta
	for(int i=0; i<fElecs.size(); ++i) {
	  fETHnt->ele[i].SetPtEtaPhiE(fTR->ElPt [fElecs[i]], fTR->ElEta[fElecs[i]], 
				      fTR->ElPhi[fElecs[i]], fTR->ElE  [fElecs[i]]); // = GetEle4Momenta(fElecs[i]);
	}
	for(int i=0; i<fMuons.size(); ++i) {
	  fETHnt->muo[i].SetPtEtaPhiM(fTR->MuPt [fMuons[i]], fTR->MuEta[fMuons[i]], 
				      fTR->MuPhi[fMuons[i]], 0.106);                     // = GetMuo4Momenta(fMuons[i]);
	}
	// Fill met
	fETHnt->pfmet[0].SetPtEtaPhiM(fTR->PFMET,0,fTR->PFMETphi,0);

	// remaining variables are filles elsewhere
	fMassTree           ->Fill();
	fETHtree            ->Fill();
}

// ************************************************************************************************
void MassAnalysis::ControlPlots(){

	// NJets
	fHNJets -> Fill(fJets.size());
	if(fIsCleanMultiJetEvent){
		fHNJets_clean -> Fill(fJets.size());
	}

	// LeptConfig
	fHLeptConfig ->Fill(fLeptConfig);
	if(fIsCleanMultiJetEvent){
		 fHLeptConfig_clean -> Fill(fLeptConfig);
	}

	
	// VSPT
	fHVectorSumPt -> Fill(fVectorSumPt);
	if(fIsCleanMultiJetEvent){
		fHVectorSumPt_clean -> Fill(fVectorSumPt); 	
	}

	// delta Phi(j,MET)
	fHDPhiJ1vsDPhiJ2 -> Fill(fDeltaPhi1, fDeltaPhi2); 	
	if(fIsCleanMultiJetEvent){
		fHDPhiJ1vsDPhiJ2_clean -> Fill(fDeltaPhi1, fDeltaPhi2); 	
	}
	// MET spectrum
	fHPFMET -> Fill(fTR->PFMET);
	if(fIsCleanMultiJetEvent){
		fHPFMET_clean -> Fill(fTR->PFMET); 
	}

	// MHT
	fHMHT -> Fill(fMHT);
	if(fIsCleanMultiJetEvent){
		fHMHT_clean -> Fill(fMHT);
	}

	// Jet Pt spectrum
	for(int i=0; i<fJets.size(); ++i){
		fHJpt  -> Fill(fTR->PFJPt[fJets[i]]);
		fHJEta -> Fill(fTR->PFJEta[fJets[i]]);
		if(fIsCleanMultiJetEvent){
			fHJpt_clean  -> Fill(fTR->PFJPt[fJets[i]]); 
			fHJEta_clean -> Fill(fTR->PFJEta[fJets[i]]);
		}
	}


	// MPT 
	fHMPT -> Fill(sqrt((fTR->TrkPtSumx)*(fTR->TrkPtSumx) + (fTR->TrkPtSumy)*(fTR->TrkPtSumy)) );
	TVector3 tracks(0.,0.,0.);
	TVector3 track(0.,0.,0.);
	for(int i=0; i< fTR->NTracks; ++i){
		track.SetPtEtaPhi(fabs(fTR->TrkPt[i]), fTR->TrkEta[i], fTR->TrkPhi[i]);
		tracks += track;
	}
	fHMPT_selected ->Fill(tracks.Pt());
		
	if(fJets.size()>1) {
		double delta_R=Util::GetDeltaR(fTR->PFJEta[0], fTR->JEta[0], fTR->PFJPhi[0], fTR->JPhi[0]);	
		fHPFandCalo_deltaR -> Fill(delta_R);
		delta_R       =Util::GetDeltaR(fTR->PFJEta[0], fTR->PFJEta[1], fTR->PFJPhi[0], fTR->PFJPhi[1]);
		fHPFJ1J2DeltaR -> Fill(delta_R);	
	}
}

// *************************************************************************************************
void MassAnalysis::PseudoJetMasses(){
	if(fJets.size() <2 ) return;

	// ---------------------------------------------
	// Cleaning

	if(fIsCleanMultiJetEvent == false) return;


	// ----------------------------------------------
	// AlphaT 
	double alphaT=-999.99;
	if( fLeptConfig==null ){
		std::vector<TLorentzVector> jets;
		for(int i=0; i<fJets.size(); ++i){
			TLorentzVector v;
	 		v.SetPxPyPzE(fTR->PFJPx[fJets[i]],fTR->PFJPy[fJets[i]], fTR->PFJPz[fJets[i]], fTR->PFJE[fJets[i]]);
	 		jets.push_back(v);
	 	}
	 	alphaT= GetAlphaT(jets);
	 	fHAlpahT_PseudoJet->Fill(alphaT);
	}
	// fill tree variable
	fTpseudojetAlphaT = alphaT;
	fETHnt->misc.PseudojetAlphaT = alphaT;




	// ---------------------------------------------
	// MT2
	//    if there are selected leptons, add them to the visible system and fill fHMT2_PseudoJetsWithLeptons
	//    if there is at least one b-jet: fill fHMT2_PseudoJetsWithB    
	//    in any case, fill fHMT2_pseudojet
	
	// make pseudojets with hemispheres
	vector<float> px, py, pz, E;
	for(int i=0; i<fJets.size(); ++i){
		px.push_back(fTR->PFJPx[fJets[i]]);
		py.push_back(fTR->PFJPy[fJets[i]]);
		pz.push_back(fTR->PFJPz[fJets[i]]);
		 E.push_back(fTR->PFJE[fJets[i]]);
	}
		
	vector<TLorentzVector> leptonmomenta = GetLepton4Momenta();
	for(int i=0; i< leptonmomenta.size(); ++i){
		px.push_back(leptonmomenta[i].Px());
		py.push_back(leptonmomenta[i].Py());
		pz.push_back(leptonmomenta[i].Pz());
		 E.push_back(leptonmomenta[i].E());
	}

	// get hemispheres (seed: max inv mass, association method: minimal lund distance)
	fHemisphere = new Hemisphere(px, py, pz, E, 2, 3);
	vector<int> grouping = fHemisphere->getGrouping();

	TLorentzVector pseudojet1(0.,0.,0.,0.);
	TLorentzVector pseudojet2(0.,0.,0.,0.);
	
	TLorentzVector pmiss;
	pmiss.SetPx(fTR->PFMETpx);
	pmiss.SetPy(fTR->PFMETpy);

	for(int i=0; i<px.size(); ++i){
		if(grouping[i]==1){
			pseudojet1.SetPx(pseudojet1.Px() + px[i]);
			pseudojet1.SetPy(pseudojet1.Py() + py[i]);
			pseudojet1.SetPz(pseudojet1.Pz() + pz[i]);
			pseudojet1.SetE( pseudojet1.E()  + E[i]);	
			fETHnt->jet[i].inHemisphere = 1;	
		}else if(grouping[i] == 2){
			pseudojet2.SetPx(pseudojet2.Px() + px[i]);
			pseudojet2.SetPy(pseudojet2.Py() + py[i]);
			pseudojet2.SetPz(pseudojet2.Pz() + pz[i]);
			pseudojet2.SetE( pseudojet2.E()  + E[i]);
			fETHnt->jet[i].inHemisphere = 2;	
		}
	}
	delete fHemisphere;
		
	// plot MT2 histos
	for(int i=0; i<fMT2_histos_number; ++i){
		double MT2       =GetMT2(pseudojet1, 0., pseudojet2, 0., pmiss, i*fMT2_histos_step);
		double MT2massive=GetMT2(pseudojet1, pseudojet1.M(), pseudojet2, pseudojet2.M(), pmiss, i*fMT2_histos_step);
		fHMT2_pseudojet[i] -> Fill(MT2);
		if(fLeptConfig!= null) {fHMT2_PseudoJetWithLeptons[i] -> Fill(MT2);}
		if(fLeptConfig== null) {fHMT2_PseudoJetNoLeptons[i]   -> Fill(MT2);}
		if(fBJets.size()>0)    {fHMT2_PseudoJetWithB[i]       -> Fill(MT2);}
	       	if(i==0)               {
			fHPseudoJetMT2vsVectorSumPt   -> Fill(MT2, fVectorSumPt);
			fHPseudoJetMT2vsMETsign       -> Fill(MT2, (fTR->PFMET)/sqrt(fTR->SumEt));
			fHPseudoJetMT2vsMET           -> Fill(MT2, fTR->PFMET);
			fHPseudoJetMT2vsAlphaT        -> Fill(MT2, alphaT);
			fHPseudoJetMT2vsLeadingJEta   -> Fill(MT2, fTR->PFJEta[fJets[0]]);
			fHPseudoJetMT2vsMHT           -> Fill(MT2, fMHT);
		}
		if(i==0 && MT2 > 240  ){
			interesting_Run.push_back(fTR->Run);
			interesting_Event.push_back(fTR->Event);
			interesting_Lumi.push_back(fTR->LumiSection);
			interesting_value.push_back(MT2);
			interesting_Type.push_back("pseudojetMT2");
			if(fSetName=="JetMET_Prompt" || fSetName=="MultiJetRunB") UserAnalysisBase::EventPrint();
		}
		if(i==0){
			// fill tree variable
			fTpseudoJetMT2        = MT2;
			fTpseudoJetMT2massive = MT2massive;
			fETHnt->misc.PseudoJetMT2 = MT2;
			if(pseudojet1.Pt() > pseudojet2.Pt()){
				fTpseudoJet1Pt = pseudojet1.Pt();
				fTpseudoJet2Pt = pseudojet2.Pt();
				fETHnt->misc.PseudoJet1Pt = pseudojet1.Pt();
				fETHnt->misc.PseudoJet2Pt = pseudojet2.Pt();
				fETHnt->pseudoJets[0].SetPxPyPzE(pseudojet1.Px(), pseudojet1.Py(),
								 pseudojet1.Pz(), pseudojet1.E ());
				fETHnt->pseudoJets[1].SetPxPyPzE(pseudojet2.Px(), pseudojet2.Py(),
								 pseudojet2.Pz(), pseudojet2.E ());
			} else {
				fTpseudoJet1Pt = pseudojet2.Pt();
				fTpseudoJet2Pt = pseudojet1.Pt();
				fETHnt->misc.PseudoJet1Pt = pseudojet2.Pt();
				fETHnt->misc.PseudoJet2Pt = pseudojet1.Pt();
				fETHnt->pseudoJets[1].SetPxPyPzE(pseudojet1.Px(), pseudojet1.Py(),
								 pseudojet1.Pz(), pseudojet1.E ());
				fETHnt->pseudoJets[0].SetPxPyPzE(pseudojet2.Px(), pseudojet2.Py(),
								 pseudojet2.Pz(), pseudojet2.E ());
			}		
		}	
	}
	// simplified MT2
	TVector3 vpseudojet1, vpseudojet2;
	vpseudojet1.SetXYZ(pseudojet1.Px(), pseudojet1.Py(), pseudojet1.Pz());
	vpseudojet2.SetXYZ(pseudojet2.Px(), pseudojet2.Py(), pseudojet2.Pz());

	fTpseudoJetMT2simple=sqrt(2*pseudojet1.Perp()*pseudojet2.Perp()*(1+cos(vpseudojet1.Angle(vpseudojet2))));


	// get the delta_phi between the two axis
	TVector3 vec1, vec2;
	vec1.SetPtEtaPhi(pseudojet1.Pt(),pseudojet1.Eta(),pseudojet1.Phi());
	vec2.SetPtEtaPhi(pseudojet2.Pt(),pseudojet2.Eta(),pseudojet2.Phi());
	fTPseudoJetMT2AxisdPhi = vec1.Angle(vec2);
	fETHnt->misc.PseudoJetMT2AxisdPhi = vec1.Angle(vec2);
	fHPseudoJetMT2AxisdPhi -> Fill(fTPseudoJetMT2AxisdPhi);	

	// ---------------------------------------------
	// MCT pseudojets
	double MCT     =GetMCT(pseudojet1, pseudojet2);
	TVector2 met;
	met.Set(pmiss.Px(), pmiss.Py());
	TLorentzVector DTM(0.,0.,0.,0.);
	double MCTcorr          =GetMCTcorr(pseudojet1, pseudojet2, DTM, met );
	TLorentzVector pseudojet1massless, pseudojet2massless;
	pseudojet1massless.SetXYZM(pseudojet1.Px(), pseudojet1.Py(), pseudojet1.Pz(), 0.);
	pseudojet2massless.SetXYZM(pseudojet2.Px(), pseudojet2.Py(), pseudojet2.Pz(), 0.);
	double MCTcorrMassless  =GetMCTcorr(pseudojet1massless, pseudojet2massless, DTM, met );
	double MCTcorrMinusMass =sqrt(MCTcorr*MCTcorr - pseudojet1.M2() -pseudojet2.M2());

	fTpseudoJetMCT          = MCTcorr;
	fTpseudoJetMCTMinusMass = MCTcorrMinusMass;
	fTpseudoJetMCTmassless  = MCTcorrMassless;
	fETHnt->misc.PseudoJetMCT  = MCTcorr;

	fHMCT_PseudoJet   -> Fill(MCT);
	if(fLeptConfig!= null) {fHMCT_PseudoJetWithLeptons     -> Fill(MCT);}
	if(fLeptConfig== null) {fHMCT_PseudoJetNoLeptons       -> Fill(MCT);}
	if(fBJets.size()>0)    {fHMCT_PseudoJetWithB           -> Fill(MCT);}
	fHMCTcorr_PseudoJet   -> Fill(MCTcorr);
	if(fLeptConfig!= null) {fHMCTcorr_PseudoJetWithLeptons -> Fill(MCTcorr);}
	if(fLeptConfig== null) {fHMCTcorr_PseudoJetNoLeptons   -> Fill(MCTcorr);}
	if(fBJets.size()>0)    {fHMCTcorr_PseudoJetWithB       -> Fill(MCTcorr);}


	if(fVerbose > 1) {
		cout << "MT2 " << fTpseudoJetMT2 << " MT2massive " << fTpseudoJetMT2massive << "  MT2simplified " << fTpseudoJetMT2simple << 
		" MCTcorr " << fTpseudoJetMCT << " MCTcorr minus masses " << fTpseudoJetMCTMinusMass << " MCTcorr massless " << MCTcorrMassless << endl;
	}



}

// **************************************************************************************************
void MassAnalysis::DiBJetMasses(){
	if(fBJets.size()!=2) return;
	
	// ---------------------------------------------
	// Cleaning
	if(fIsCleanJetEvent == false) return;

	TLorentzVector p1, p2;
	p1.SetPxPyPzE(fTR->JPx[fBJets[0]], fTR->JPy[fBJets[0]], fTR->JPz[fBJets[0]], fTR->JE[fBJets[0]]);
	p2.SetPxPyPzE(fTR->JPx[fBJets[1]], fTR->JPy[fBJets[1]], fTR->JPz[fBJets[1]], fTR->JE[fBJets[1]]);

	TLorentzVector pmiss(0., 0., 0., 0.);
	pmiss.SetPx(fTR->PFMETpx);
	pmiss.SetPy(fTR->PFMETpy);

	// MT2 in case of additional OS di-leptons
	// di-b-jets with leptons summed up to MET
  	// for TTbar events, the parent mass would then be the top!
	if(fLeptConfig == OS_ee || fLeptConfig == OS_emu ||fLeptConfig == OS_mumu ){
		vector<TLorentzVector> leptmomenta = GetLepton4Momenta();
		TLorentzVector p_invisible=pmiss + leptmomenta[0] + leptmomenta[1];
		for(int i=0; i<fMT2_histos_number; ++i){	
			double MT2=GetMT2(p1, 0., p2, 0., p_invisible, i*fMT2_histos_step);
			fHMT2_diBjetdiLept[i] -> Fill(MT2);
		}
	}

	// -------------------------------------------
	// lepton veto
	if(fLeptConfig !=null) return;

	// -------------------------------------------
	// invariant mass
	double mass = (p1+p2).M();
	fHInvMassDiBjet -> Fill(mass);

	// -------------------------------------------
	// MCT
	double MCT=GetMCT(p1, p2);
	fHMCT_DiBJet -> Fill(MCT);
	// MCTcorr
	TVector2 met(fTR->PFMETpx, fTR->PFMETpy);
	TLorentzVector DTM(0.,0.,0.,0.);
	double MCTcorr = GetMCTcorr(p1, p2, DTM, met);
	fHMCTcorr_DiBJet -> Fill(MCTcorr);

	// ------------------------------------------
	// MT2
	for(int i=0; i<fMT2_histos_number; ++i){
		double MT2=GetMT2(p1, 0., p2, 0., pmiss, i*fMT2_histos_step);
		fHMT2_diBjet[i] ->Fill(MT2);
	}
	
}

// ***************************************************************************************************
void MassAnalysis::MultiBMasses(){
	if( fBJets.size()<3 ) return;

	// inv mass of two bjets if they are in the same hemisphere
	// the idea: in SUSY, higgs from decay of heavy particle -> boosted. H->b bbar, with two b in same hemi
	vector<float> px, py, pz, E;
	for(int i=0; i<fJets.size(); ++i ){
		px.push_back(fTR->JPx[fJets[i]]);
		py.push_back(fTR->JPy[fJets[i]]);
		pz.push_back(fTR->JPz[fJets[i]]);
		E.push_back(fTR->JE[fJets[i]]);
	}
		
	fHemisphere = new Hemisphere(px, py, pz, E, 2, 3);
	vector<int> grouping = fHemisphere->getGrouping();
	vector<int> b_group1, b_group2;
	for(int i=0; i<fJets.size(); ++i){
		if(grouping[i]==1 && fTR->JbTagProbSimpSVHighPur[fJets[i]] > 2.0) b_group1.push_back(fJets[i]);
		if(grouping[i]==2 && fTR->JbTagProbSimpSVHighPur[fJets[i]] > 2.0) b_group2.push_back(fJets[i]);
	}

	TLorentzVector b1, b2;
	if(b_group1.size()==2 && b_group2.size()<2){
		b1.SetPxPyPzE(fTR->JPx[b_group1[0]],fTR->JPy[b_group1[0]],fTR->JPz[b_group1[0]],fTR->JE[b_group1[0]]);
		b2.SetPxPyPzE(fTR->JPx[b_group1[1]],fTR->JPy[b_group1[1]],fTR->JPz[b_group1[1]],fTR->JE[b_group1[1]]);
	} else if(b_group1.size()<2 && b_group2.size()==2){
		b1.SetPxPyPzE(fTR->JPx[b_group2[0]],fTR->JPy[b_group2[0]],fTR->JPz[b_group2[0]],fTR->JE[b_group2[0]]);
		b2.SetPxPyPzE(fTR->JPx[b_group2[1]],fTR->JPy[b_group2[1]],fTR->JPz[b_group2[1]],fTR->JE[b_group2[1]]);
	}

	double invmass = (b1+b2).M();
		
	fHInvMassdiBHemi -> Fill(invmass);

}



// ***************************************************************************************************
void MassAnalysis::DiJetMasses(){
	if( fJets.size()!=2 ) return;    // only look at dijets
	if( fIsCleanJetEvent == false) return;  // reject events that fails the cleaning requirements
	
	// -------------------------------------------
	// lepton veto
	if(fLeptConfig !=null) return;

	// filling of jets 
	TLorentzVector p1;
	TLorentzVector p2;		
	p1.SetPxPyPzE(fTR->PFJPx[fJets[0]], fTR->PFJPy[fJets[0]], fTR->PFJPz[fJets[0]], fTR->PFJE[fJets[0]]);
	p2.SetPxPyPzE(fTR->PFJPx[fJets[1]], fTR->PFJPy[fJets[1]], fTR->PFJPz[fJets[1]], fTR->PFJE[fJets[1]]);
	
	// --------------------------------------------------
	// invariant mass
	fHInvMassDijet -> Fill((p1+p2).M());

	// ------------------------------------------------
	// MCT
	double MCT=GetMCT(p1, p2);
	fHMCT_DiJet -> Fill(MCT);

	// MCTcorr
	TVector2 met(fTR->PFMETpx, fTR->PFMETpy);
	TLorentzVector DTM(0.,0.,0.,0.);
	double MCTcorr = GetMCTcorr(p1, p2, DTM, met);
	fHMCTcorr_DiJet ->Fill(MCTcorr);

	if(MCT > 240){
		interesting_Run.push_back(fTR->Run);
		interesting_Event.push_back(fTR->Event);
		interesting_Lumi.push_back(fTR->LumiSection);
		interesting_value.push_back(MCT);
		interesting_Type.push_back("dijetMCT");
		if(fTR->Run == 143657 && fTR->LumiSection == 1310 && fTR->Event==1285062562) {
			PrintEvent();
			cout << "dijet mass " << (p1+p2).M() << endl;
		}
	}

	// -------------------------------------------------
	// MT2
	TLorentzVector pmiss(0., 0., 0., 0.);
	pmiss.SetPx(fTR->PFMETpx);
	pmiss.SetPy(fTR->PFMETpy);
	for(int i=0; i<fMT2_histos_number; ++i){
		double MT2=GetMT2(p1, 0., p2, 0., pmiss, i*fMT2_histos_step);
		fHMT2_dijet[i] ->Fill(MT2);
		if(i==0){fHVectorSumPtvsDiJetMT2->Fill(fVectorSumPt, MT2);}
		if(i==0 && MT2 > 240){
			interesting_Run.push_back(fTR->Run);
			interesting_Event.push_back(fTR->Event);
			interesting_Lumi.push_back(fTR->LumiSection);
			interesting_value.push_back(MT2);
			interesting_Type.push_back("dijetMT2");

		}
	}

	// --------------------------------------------------
	// alphaT  
	if( fTR->PFJPt[fJets[0]] > 150 ){
		std::vector<TLorentzVector> jets;
		jets.push_back(p1);
		jets.push_back(p2);
	 	double alphaT= GetAlphaT(jets);
	 	fHAlpahT_DiJet->Fill(alphaT);
	}
}

// ****************************************************************************************************
void MassAnalysis::SingleLeptonMasses(){
	if(! (fLeptConfig==e || fLeptConfig==mu) ) return;
	if(fTR->PFMET < 30) return;

	// compute transverse mass
	TLorentzVector l;
	if(fLeptConfig==e) {l.SetPtEtaPhiE(fTR->ElPt[fElecs[0]],fTR->ElEta[fElecs[0]],fTR->ElPhi[fElecs[0]],fTR->ElE[fElecs[0]]);}
	if(fLeptConfig==mu){l.SetPtEtaPhiM(fTR->MuPt[fMuons[0]],fTR->MuEta[fMuons[0]],fTR->MuPhi[fMuons[0]],0.105);}

	double ETlept = sqrt(l.M2() + l.Perp2());
	double METpt  = fTR->PFMET;
	double METpx  = fTR->PFMETpx;
	double METpy  = fTR->PFMETpy;
	double MT     = sqrt(2*METpt*ETlept - 2*l.Px()*METpx - 2*l.Py()*METpy );

	if(fLeptConfig==e ){fHMT_single_e -> Fill(MT);}
	if(fLeptConfig==mu){fHMT_single_mu-> Fill(MT);}

	
	if(fJets.size()!=0) return; // plots also MT in case there are no selected jets (pure W -> l nu selection..)

	if(fLeptConfig==e ){fHMT_single_e_nojets -> Fill(MT);}
	if(fLeptConfig==mu){fHMT_single_mu_nojets-> Fill(MT);}

}

// ****************************************************************************************************
void MassAnalysis::DiLeptonMasses(){
	if(! (fLeptConfig==OS_ee || fLeptConfig==OS_mumu || fLeptConfig==OS_emu || fLeptConfig==SS_ee || fLeptConfig==SS_mumu || fLeptConfig==SS_emu) ) return;

	// TLorentzVector of di-leptons
	vector<TLorentzVector> dilept4mom = GetLepton4Momenta();
	TLorentzVector p1=dilept4mom[0];
	TLorentzVector p2=dilept4mom[1];

	// define Upstream Transverse Momentum
	TLorentzVector P_UTM(0.,0.,0.,0.);
	for(int i=0; i<fJets.size(); ++i){
		TLorentzVector p;
		p.SetPxPyPzE(fTR->PFJPx[fJets[i]],fTR->PFJPy[fJets[i]],fTR->PFJPz[fJets[i]],fTR->PFJE[fJets[i]]);
		P_UTM +=p;
	}

	// calculate inv mass, MCT, MCTperp, MT2, MT2perp	
	double MCTperp=-1.;
	double MCTcorr=-1.;
	double MT2perp[fMT2_histos_number];
	double MT2[fMT2_histos_number];

	double MCT     =GetMCT(p1, p2);
	double mass    =(p1+p2).M();

	// check Zmass window for mass plots other than invmass
	bool Zmassveto(false);
	if(fLeptConfig==OS_ee ||fLeptConfig==OS_mumu ){
		if( (mass > fCut_DiLeptOSSFInvMass_lowercut) && (mass < fCut_DiLeptOSSFInvMass_uppercut) ) {Zmassveto=true;}
	}

	
	// calculate mass variables
	if(P_UTM != (0.,0.,0.,0.)) {
		MCTperp =GetMCTperp(p1, p2, P_UTM);
		MCTcorr =GetMCTcorr(p1, p2, P_UTM);
	}
	TLorentzVector pmiss(0., 0., 0., 0.);
	pmiss.SetPx(fTR->PFMETpx);
	pmiss.SetPy(fTR->PFMETpy);
	for(int i=0; i<fMT2_histos_number; ++i){
		int m_invisible = i*fMT2_histos_step;
		MT2[i]    =GetMT2(p1, p1.M(), p2, p2.M(), pmiss, m_invisible);
		if(P_UTM != (0.,0.,0.,0.)){
			MT2perp[i]=GetMT2perp(p1, p2, P_UTM, m_invisible);
		}else   MT2perp[i]=-1.;
	}
		
		
	// fill histos
	if(fLeptConfig==OS_ee){
		fHInvMassOSee   -> Fill(mass);
		fHInvMassOSll   -> Fill(mass);	
		if(Zmassveto==false){
			if(P_UTM != (0.,0.,0.,0.)) {
				fHMCTperp_OSee -> Fill(MCTperp);
				fHMCTcorr_OSee -> Fill(MCTcorr);
			}
			fHMCT_OSee      -> Fill(MCT);
			for(int i=0; i<fMT2_histos_number; ++i){
				fHMT2_OSll[i]     -> Fill(MT2[i]);
				fHMT2_OSee[i]     -> Fill(MT2[i]);
				if(P_UTM != (0.,0.,0.,0.)) {
					fHMT2perp_OSee[i] -> Fill(MT2perp[i]);
				}	
			}
		}
	} else if(fLeptConfig==OS_mumu){
		fHInvMassOSmumu -> Fill(mass);
		fHInvMassOSll   -> Fill(mass);	
		if(Zmassveto==false){
			if(P_UTM != (0.,0.,0.,0.)) {
				fHMCTperp_OSmumu -> Fill(MCTperp);
				fHMCTcorr_OSmumu -> Fill(MCTcorr);
			}
			fHMCT_OSmumu      -> Fill(MCT);
			for(int i=0; i<fMT2_histos_number; ++i){
				fHMT2_OSll[i]       -> Fill(MT2[i]);
				fHMT2_OSmumu[i]     -> Fill(MT2[i]);
				if(P_UTM != (0.,0.,0.,0.)) {
					fHMT2perp_OSmumu[i] -> Fill(MT2perp[i]);
				}
			}
		}
	} else if(fLeptConfig==OS_emu){
		fHInvMassOSemu  -> Fill(mass);
		fHInvMassOSll   -> Fill(mass);	
		if(P_UTM != (0.,0.,0.,0.)) {
			fHMCTperp_OSemu -> Fill(MCTperp);
			fHMCTcorr_OSemu -> Fill(MCTcorr);
		}
		fHMCT_OSemu     -> Fill(MCT);
		for(int i=0; i<fMT2_histos_number; ++i){
			fHMT2_OSll[i]       -> Fill(MT2[i]);
			fHMT2_OSemu[i]      -> Fill(MT2[i]);
			if(P_UTM != (0.,0.,0.,0.)) {
				fHMT2perp_OSemu[i]  -> Fill(MT2perp[i]);
			}
		}
	} else if(fLeptConfig==SS_emu || fLeptConfig==SS_ee || fLeptConfig==SS_mumu ){
		for(int i=0; i<fMT2_histos_number; ++i){
			fHMT2_SSll[i]    -> Fill(MT2[i]);
		}
		fHInvMassSSll -> Fill(mass);
		if(fLeptConfig==SS_emu ) fHInvMassSSemu  -> Fill(mass);
		if(fLeptConfig==SS_ee  ) fHInvMassSSee   -> Fill(mass);
		if(fLeptConfig==SS_mumu) fHInvMassSSmumu -> Fill(mass);
	}
}

 // ****************************************************************************************************
double MassAnalysis::GetMCTperp(TLorentzVector p1, TLorentzVector p2, TLorentzVector P_UTM){
	double m1=p1.M();
	double m2=p2.M();
	
	TVector3 p1_T_perp  =  GetMomPerp(p1, P_UTM); 
	TVector3 p2_T_perp  =  GetMomPerp(p2, P_UTM);
	
	double E1_T_perp = sqrt( pow(m1,2) + p1_T_perp.Mag2());
	double E2_T_perp = sqrt( pow(m2,2) + p2_T_perp.Mag2());	
	
//	double MCT_perp= sqrt( pow(m1,2) + pow(m2,2) + 2 *(E1_T_perp*E2_T_perp + p1_T_perp.Dot(p2_T_perp)) );
	double MCT_perp= sqrt( pow(m1,2) + pow(m2,2) + 2 *(E1_T_perp*E2_T_perp + p1_T_perp.Mag()*p2_T_perp.Mag()*cos(p1_T_perp.Angle(p2_T_perp))) );


	// compare to Tovey Impelemntation
	TLorentzVector DTM(0.,0.,0.,0.);
	TVector2       pmiss;
	pmiss.Set(-P_UTM.Px()-p1.Px()-p2.Px(), -P_UTM.Py()-p1.Py()-p2.Py());
	double ToveyMCTperp=GetToveyMCTperp(p1, p2, DTM, pmiss);
	// if(fabs(MCT_perp - ToveyMCTperp)/ToveyMCTperp > 0.01) cout << "MCTperp" << MCT_perp << " Tovey " << ToveyMCTperp << endl;


	return MCT_perp;
}

// *****************************************************************************************************
double MassAnalysis::GetToveyMCTperp(TLorentzVector p1, TLorentzVector p2, TLorentzVector DTM, TVector2 pmiss){
	fMCT = new TMctLib;
	double MCTperp = fMCT -> mcy(p1, p2, DTM, pmiss);
	delete fMCT;
	return MCTperp;
}

// *****************************************************************************************************
double MassAnalysis::GetMCTcorr(TLorentzVector p1, TLorentzVector p2, TLorentzVector DTM, TVector2 pmiss){
	fMCT = new TMctLib;
	double MCTcorr = fMCT -> mctcorr(p1, p2, DTM, pmiss, 7000., 0.);
	delete fMCT;
	return MCTcorr;
}

// ******************************************************************************************************
double MassAnalysis::GetMCTcorr(TLorentzVector p1, TLorentzVector p2, TLorentzVector UTM){
	TVector2 pmiss;
	pmiss.Set(-UTM.Px()-p1.Px()-p2.Px(), -UTM.Py()-p1.Py()-p2.Py());
	TLorentzVector DTM(0.,0.,0.,0.);

	fMCT = new TMctLib;
	double MCTcorr = fMCT -> mctcorr(p1, p2, DTM, pmiss, 7000., 0.);
	delete fMCT;
	return MCTcorr;

}
// ******************************************************************************************************
double MassAnalysis::GetMT2perp(TLorentzVector p1, TLorentzVector p2, TLorentzVector P_UTM, double m_inv){
	TVector3 p1_T_perp =  GetMomPerp(p1, P_UTM);
	TVector3 p2_T_perp =  GetMomPerp(p2, P_UTM);
	
	// double A_T_perp   = 1/2.*(p1_T_perp.Mag()*p2_T_perp.Mag()+p1_T_perp.Dot(p2_T_perp) );
	// above implementation is problematic: A_T_perp is sometimes <0 (presicion problem) 
	// workaround: implementation below
	double A_T_perp   = 1/2.*( p1_T_perp.Mag()*p2_T_perp.Mag() * (1+cos(p1_T_perp.Angle(p2_T_perp)) ));
	double MT2_perp   = sqrt(A_T_perp) + sqrt(A_T_perp + pow(m_inv,2));
	
	return MT2_perp;
}

// *******************************************************************************************************
TVector3 MassAnalysis::GetMomPerp(TLorentzVector p, TLorentzVector P_UTM){
	TVector3 p_T(p.Px(), p.Py(), 0.);
	TVector3 pUTM_T(P_UTM.Px(), P_UTM.Py(), 0.);
	double pUTM_T_norm2=pUTM_T.Mag2();
	
	TVector3 p_T_perp = p_T - 1/pUTM_T_norm2*(p_T.Dot(pUTM_T))*pUTM_T ;	
	return p_T_perp;
}

// ****************************************************************************************************
void MassAnalysis::MassesForTTbar(){
	// MT2 and MCTperp for TTbar events using Gen Info
	vector<TLorentzVector> leptons;
	for(int i=0; i<fTR->NGenLeptons; ++i){
		if(! (abs(fTR->GenLeptonID[i]) == 11 || abs(fTR->GenLeptonID[i]) == 13) ) continue;
		if(! (abs(fTR->GenLeptonMID[i]) == 24) ) continue;
		if(! (abs(fTR->GenLeptonGMID[i]) == 6) ) continue;
		TLorentzVector p;
		p.SetPtEtaPhiM(fTR->GenLeptonPt[i],fTR->GenLeptonEta[i],fTR->GenLeptonPhi[i],0.);
		leptons.push_back(p);
	}
	vector<TLorentzVector> neutrinos;
	for(int i=0; i<fTR->NGenLeptons; ++i){
		if(! (abs(fTR->GenLeptonID[i]) == 12 || abs(fTR->GenLeptonID[i]) == 14) ) continue;
		if(! (abs(fTR->GenLeptonMID[i]) == 24) ) continue;
		if(! (abs(fTR->GenLeptonGMID[i]) == 6) ) continue;
		TLorentzVector p;
		p.SetPtEtaPhiM(fTR->GenLeptonPt[i],fTR->GenLeptonEta[i],fTR->GenLeptonPhi[i],0.);
		neutrinos.push_back(p);
	}
	if(leptons.size()==2 && neutrinos.size()==2){
		// define P_UTM: negative vector sum of leptons and neutrinos
		//               this should be equal to the two b-jets (no quarks in genparticles)
		double Px_UTM = -(leptons[0].Px()+leptons[1].Px()+neutrinos[0].Px()+neutrinos[1].Px() );
		double Py_UTM = -(leptons[0].Py()+leptons[1].Py()+neutrinos[0].Py()+neutrinos[1].Py() );
		TLorentzVector P_UTM(0., 0., 0., 0.);
		P_UTM.SetPx(Px_UTM);
		P_UTM.SetPy(Py_UTM);
		
		// MCTperp, MT2perp, MCT, MT2
		TLorentzVector pmiss = neutrinos[0]+neutrinos[1];
		double MT2    =GetMT2(leptons[0], leptons[0].M(), leptons[1], leptons[1].M(), pmiss, 0.);
		double MCT    =GetMCT(leptons[0], leptons[1]);
		double MCTperp=GetMCTperp(leptons[0], leptons[1], P_UTM);
		double MT2perp=GetMT2perp(leptons[0], leptons[1], P_UTM, 0.);
		double MCTcorr=GetMCTcorr(leptons[0], leptons[1], P_UTM);
		fHMCT_TTbar     -> Fill(MCT);
		fHMCTperp_TTbar -> Fill(MCTperp);
		fHMT2perp_TTbar -> Fill(MT2perp);
		fHMT2_TTbar     -> Fill(MT2);
		fHMCTcorr_TTbar -> Fill(MCTcorr);
	}
}

// ****************************************************************************************************
double MassAnalysis::GetMCT(TLorentzVector p1, TLorentzVector p2){
	double m1=p1.M();
	double m2=p2.M();
	double ET1=sqrt( pow(p1.Pt(),2) +pow(m1,2) );
	double ET2=sqrt( pow(p2.Pt(),2) +pow(m2,2) );
	
	double MCTsquared=pow(m1,2) + pow(m2,2) +2*(ET1*ET2+ (p1.Px() * p2.Px()) + (p1.Py() * p2.Py()));
	
	// Tovey implementation
	fMCT = new TMctLib();
	double Tovey_MCT = fMCT -> mct(p1, p2);
	delete fMCT;
//	if(fabs(Tovey_MCT -sqrt(MCTsquared))/Tovey_MCT > 0.01) cout << "Tovey_MCT "<< Tovey_MCT << " MCT " << sqrt(MCTsquared) << endl;

	return sqrt(MCTsquared);
}

// ****************************************************************************************************
double MassAnalysis::GetAlphaT(std::vector<TLorentzVector>& p4s) {
    if(p4s.size()<2) return 0;
    
    std::vector<double> pTs; for(unsigned i=0; i<p4s.size(); i++) pTs.push_back(p4s[i].Pt());
    const std::vector<double> DHTper( DeltaSumPt_permutations(p4s) );
    
    const double mDHT = *(std::min_element( DHTper.begin(), DHTper.end() ));
    const double sumPT = accumulate( pTs.begin(), pTs.end(), double(0) );
    const TLorentzVector sumP4 = accumulate( p4s.begin(), p4s.end(), TLorentzVector() );

    return 0.5 * ( sumPT - mDHT ) / sqrt( sumPT*sumPT - sumP4.Perp2() );
}

std::vector<double> MassAnalysis::DeltaSumPt_permutations(std::vector<TLorentzVector>& p4s) {
  	std::vector<std::vector<double> > ht( 1<<(p4s.size()-1) , std::vector<double>( 2, 0.) ); // initializiert einen vector ht der size 1<<(p4s.size()-1), wobei jeder eintrag ein vector (0,0) ist
  	for(unsigned i=0; i < ht.size(); i++) {
    	for(unsigned j=0; j < p4s.size(); j++) {
			ht [i] [(i/(1<<j))%2] += p4s[j].Pt();
    	}
  	}
  	std::vector<double> deltaHT; for(unsigned i=0; i<ht.size(); i++) deltaHT.push_back(fabs(ht[i][0]-ht[i][1]));
  	return deltaHT;
}


// ****************************************************************************************************
double MassAnalysis::GetMT2(TLorentzVector v1, double mv1, TLorentzVector v2, double mv2, TLorentzVector p_unobs, int m_invisible){
	double pa[3];
	double pb[3];
	double pmiss[3];

	pmiss[1]=p_unobs.Px();
	pmiss[2]=p_unobs.Py();

	pa[0]=mv1;
	pa[1]=v1.Px();
	pa[2]=v1.Py();

	pb[0]=mv2;
	pb[1]=v2.Px();
	pb[2]=v2.Py();
	
	fMT2 = new Davismt2();
	fMT2->set_momenta(pa, pb, pmiss);
	fMT2->set_mn(m_invisible);
	double MT2=fMT2->get_mt2();
	delete fMT2;
	return MT2;
}


// ****************************************************************************************************
vector<TLorentzVector> MassAnalysis::GetLepton4Momenta(){
	vector<TLorentzVector> momenta;
	for(int i=0; i<fElecs.size(); ++i){
		TLorentzVector p;
		p.SetPtEtaPhiE(fTR->ElPt[fElecs[i]],fTR->ElEta[fElecs[i]], fTR->ElPhi[fElecs[i]], fTR->ElE[fElecs[i]]);
		momenta.push_back(p);
	}
	for(int i=0; i<fMuons.size(); ++i){
		TLorentzVector p;
		p.SetPtEtaPhiM(fTR->MuPt[fMuons[i]],fTR->MuEta[fMuons[i]], fTR->MuPhi[fMuons[i]],0.105);
		momenta.push_back(p);
	}
	return momenta;
}


// ***************************************************************************************************
void MassAnalysis::SetWeight(float weight){
	fWeight = weight;
}


// ***************************************************************************************************
void MassAnalysis::PrintEvent(){
	cout << "********************************************************************************" << endl;
	cout << "PrintEvent for " << fTR->Run << ":" << fTR->LumiSection << ":" << fTR->Event      << endl;
	cout << endl;
	for(int i=0; i<fJets.size(); ++i){
		cout << "Jet" << i+1 << ": Pt " << fTR->PFJPt[fJets[i]]  
		     << " Phi " << fTR->PFJPhi[fJets[i]] << " Eta " << fTR->PFJEta[fJets[i]] <<endl;
	}
	for(int i=0; i<fElecs.size(); ++i){
		cout << "El"  << i+1 << ": Pt " << fTR->ElPt[fElecs[i]]
		     << " Phi " << fTR->ElPhi[fElecs[i]] << " Eta " << fTR->ElPhi[fElecs[i]] << endl;
	}
	for(int i=0; i<fMuons.size(); ++i){
		cout << "Mu"  << i+1 << ": Pt " << fTR->MuPt[fMuons[i]]
		     << " Phi " << fTR->MuPhi[fMuons[i]] << " Eta " << fTR->MuPhi[fMuons[i]] << endl;
	}
	cout << "ptvectorsum " << fVectorSumPt << endl;
	cout << "PFMET        " << fTR->PFMET          << " PFMETphi         " << fTR->PFMETphi        << endl;
	cout << "MuJESCorrMET " << fTR->MuJESCorrMET   << " MuJESCorrMETphi  " << fTR->MuJESCorrMETphi << endl;
	cout << "DeltaPhi1    " << fDeltaPhi1   << " DeltaPhi2 " << fDeltaPhi2 << endl;
	cout << "R12          " << fR12  << " fR21 " << fR21 << endl;
	cout << "*********************************************************************************" << endl;
}
// ****************************************************************************************************
void MassAnalysis::End(){
	// cout interesting events
	cout << " *************************************************************** " << endl;
	cout << " interesting events                                              " << endl;
	for(int i=0; i< interesting_Event.size(); ++i){
		cout << interesting_Type[i] << "Run:Lumi:Evt" << " " << interesting_Run[i] << ":" << interesting_Lumi[i]
		     << ":" << interesting_Event[i] << " value " << interesting_value[i] << endl;	     
	}
	cout << " *************************************************************** " << endl;

	
	fHistFile->cd();	

	// write tree
	fMassTree->Write();
	fETHtree->Write();
		
	// fill OS ll minus emu MT2 histo
	for(int i=0; i<fMT2_histos_number; ++i){
		for(int bin=1; bin<= fHMT2_OSee[0]->GetNbinsX(); ++bin){
			double bincont = fHMT2_OSee[i]->GetBinContent(bin)+fHMT2_OSmumu[i]->GetBinContent(bin);
			fHMT2_OSllminusemu[i] ->SetBinContent(bin, bincont);
		}
	}
	for(int i=0; i<fMT2_histos_number; ++i){
		fHMT2_OSllminusemu[i]->Write();
		fHMT2_OSll[i]        ->Write();
		fHMT2_OSee[i]        ->Write();
		fHMT2_OSmumu[i]      ->Write();
		fHMT2_OSemu[i]       ->Write();
		
		fHMT2_SSll[i]        ->Write();
	
		fHMT2_pseudojet[i]            ->Write();
		fHMT2_PseudoJetWithB[i]       ->Write();
		fHMT2_PseudoJetWithLeptons[i] ->Write();
		fHMT2_PseudoJetNoLeptons[i]    ->Write();
	        
		fHMT2_diBjetdiLept[i]->Write();
		fHMT2_dijet[i]       ->Write();
		fHMT2_diBjet[i]      ->Write();
		
		fHMT2perp_OSee[i]    ->Write();
		fHMT2perp_OSmumu[i]  ->Write();
		fHMT2perp_OSemu[i]   ->Write();
	}
	fHMT2_vs_M_OSDiLept      ->Write();
	
	fHAlpahT_DiJet           ->Write();
	fHAlpahT_PseudoJet       ->Write();
 
	fHMCT_PseudoJet          ->Write();
	fHMCT_PseudoJetWithLeptons ->Write();
	fHMCT_PseudoJetNoLeptons ->Write();
	fHMCT_PseudoJetWithB     ->Write();
	fHMCT_DiJet              ->Write();
	fHMCT_DiBJet             ->Write();
	fHMCT_OSee               ->Write();
	fHMCT_OSmumu             ->Write();
	fHMCT_OSemu              ->Write();
	
	fHMCTcorr_PseudoJet          ->Write();
	fHMCTcorr_PseudoJetWithLeptons ->Write();
	fHMCTcorr_PseudoJetNoLeptons ->Write();
	fHMCTcorr_PseudoJetWithB     ->Write();
	fHMCTcorr_DiJet              ->Write();
	fHMCTcorr_DiBJet             ->Write();
	fHMCTcorr_OSee               ->Write();
	fHMCTcorr_OSmumu             ->Write();
	fHMCTcorr_OSemu              ->Write();


	fHMCTperp_OSee           ->Write();
	fHMCTperp_OSmumu         ->Write();
	fHMCTperp_OSemu          ->Write();

	fHMCTcorr_TTbar          ->Write();	
	fHMCTperp_TTbar          ->Write();
	fHMT2_TTbar              ->Write();
	fHMT2perp_TTbar          ->Write();
	fHMCT_TTbar              ->Write();
	
	fHInvMassOSee            ->Write();
	fHInvMassOSmumu          ->Write();
	fHInvMassOSemu           ->Write();
	fHInvMassOSll            ->Write();
	fHInvMassSSee            ->Write();
	fHInvMassSSmumu          ->Write();
	fHInvMassSSemu           ->Write();
	fHInvMassSSll            ->Write();
	fHInvMassDijet           ->Write();
	fHInvMassDiBjet          ->Write(); 
	fHInvMassdiBHemi         ->Write();

	fHDPhiJ1vsDPhiJ2         ->Write();
	fHDPhiJ1vsDPhiJ2_clean   ->Write();

	fHMPT                    ->Write();
	fHMPT_selected           ->Write();
	fHMHT                    ->Write();
	fHMHT_clean              ->Write();
	fHPFMET                  ->Write();	
	fHPFMET_clean            ->Write();	
	fHJpt_clean              ->Write();  
	fHJpt                    ->Write();  
	fHJEta_clean             ->Write();
	fHJEta                   ->Write();
	fHBJpt                   ->Write();
	fHBJEta                  ->Write();
	fHNJets                  ->Write();
	fHNJets_clean            ->Write();
	fHLeptConfig             ->Write();
	fHLeptConfig_clean       ->Write();
	fHVectorSumPt_clean      ->Write();
	fHVectorSumPt            ->Write();
	fHPseudoJetMT2AxisdPhi   ->Write();

	fHVectorSumPtvsDiJetMT2  ->Write();


	fHPseudoJetMT2vsVectorSumPt  ->Write();
	fHPseudoJetMT2vsMETsign      ->Write(); 
	fHPseudoJetMT2vsMET          ->Write(); 
	fHPseudoJetMT2vsAlphaT       ->Write(); 
	fHPseudoJetMT2vsLeadingJEta  ->Write();
       	fHPseudoJetMT2vsMHT	     ->Write();


	fHMT_single_e            ->Write();
	fHMT_single_mu           ->Write(); 
	fHMT_single_e_nojets     ->Write();
	fHMT_single_mu_nojets    ->Write();

	fHPFandCalo_deltaR       ->Write();
	fHPFJ1J2DeltaR           ->Write();
	fHistFile                ->Close();
}
