#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/JPTJetCollection.h"
#include "DataFormats/JetReco/interface/JetID.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "DataFormats/JetReco/interface/JetTracksAssociation.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "DiLeptonAnalysis/NTupleProducer/interface/JetFillerBase.h"


//________________________________________________________________________________________
JetFillerBase::JetFillerBase( const edm::ParameterSet& cfg, TTree* tree, 
                              const bool& isRealData )
  : FillerBase(cfg,tree,isRealData)
{
	
  // Retrieve configuration parameters
  std::string jettype = cfg.getUntrackedParameter<edm::InputTag>("tag").label();

  // parse fTag label for jet-type
  if      (std::string::npos != jettype.find("PF"))   setJetType(PF);
  else if (std::string::npos != jettype.find("Calo")) setJetType(CALO);
  else if (std::string::npos != jettype.find("JPT"))  setJetType(JPT);  
  else {
    setJetType(unknown);
    edm::LogWarning("NTP") << "!! Don't know JetType !!" << jettype;
  }


  //FIXME: could be set from configuration file...
  gMaxnobjs = 100;

  // Define all arrays
  fTpx                     = new double[gMaxnobjs];
  fTpy                     = new double[gMaxnobjs];
  fTpz                     = new double[gMaxnobjs];
  fTpt                     = new double[gMaxnobjs];
  fTe                      = new double[gMaxnobjs];
  fTet                     = new double[gMaxnobjs];
  fTeta                    = new double[gMaxnobjs];
  fTphi                    = new double[gMaxnobjs];
  fTscale                  = new double[gMaxnobjs];
  fTL1FastJetScale         = new double[gMaxnobjs];
  fTarea                   = new double[gMaxnobjs];
  fTjbTagProbTkCntHighEff  = new double[gMaxnobjs];
  fTjbTagProbTkCntHighPur  = new double[gMaxnobjs];
  fTjbTagProbSimpSVHighEff = new double[gMaxnobjs];
  fTjbTagProbSimpSVHighPur = new double[gMaxnobjs];
  fTIDLoose               = new int[gMaxnobjs]; 
  if(fJetType==CALO) {
    fTID_HPD        = new double[gMaxnobjs];
    fTID_RBX        = new double[gMaxnobjs];
    fTID_n90Hits    = new double[gMaxnobjs];
    fTID_resEMF     = new double[gMaxnobjs];
    fTEMfrac        = new double[gMaxnobjs];
    fTjChfrac       = new double[gMaxnobjs];
    fTjnAssoTracks  = new int[gMaxnobjs];
    fTNConstituents = new int[gMaxnobjs];		
    fTn90           = new int[gMaxnobjs];
  } else if (fJetType==JPT) {
    fTChMult        = new int[gMaxnobjs];
    fTID_HPD        = new double[gMaxnobjs];
    fTID_RBX        = new double[gMaxnobjs];
    fTID_n90Hits    = new double[gMaxnobjs];
    fTID_resEMF     = new double[gMaxnobjs];
  } else if (fJetType==PF) {
    fTChMult        = new int[gMaxnobjs];
    fTChHadFrac     = new double[gMaxnobjs];
    fTNeuHadFrac    = new double[gMaxnobjs];
    fTChEmFrac      = new double[gMaxnobjs];
    fTNeuEmFrac     = new double[gMaxnobjs];
    fTChMuFrac      = new double[gMaxnobjs];
    fTNeuMult       = new int[gMaxnobjs];	
    fTNConstituents = new int[gMaxnobjs];
  }
}


//________________________________________________________________________________________
JetFillerBase::~JetFillerBase(void) {
  
  // Delete all arrays
  delete [] fTpx;
  delete [] fTpy;
  delete [] fTpz;
  delete [] fTpt;
  delete [] fTe;
  delete [] fTet;
  delete [] fTeta;
  delete [] fTarea;
  delete [] fTphi;
  delete [] fTscale;
  delete [] fTL1FastJetScale;
  delete [] fTjbTagProbTkCntHighEff ;
  delete [] fTjbTagProbTkCntHighPur ;
  delete [] fTjbTagProbSimpSVHighEff;
  delete [] fTjbTagProbSimpSVHighPur;
  delete [] fTIDLoose; 
  if (fJetType==CALO) {
    delete [] fTID_HPD;
    delete [] fTID_RBX; 
    delete [] fTID_n90Hits;
    delete [] fTID_resEMF;
    delete [] fTEMfrac;
    delete [] fTjnAssoTracks;
    delete [] fTjChfrac;
    delete [] fTn90;		
    delete [] fTNConstituents;		
  } else if (fJetType==JPT) {
    delete [] fTChMult;     
    delete [] fTID_HPD;    
    delete [] fTID_RBX;    
    delete [] fTID_n90Hits; 
    delete [] fTID_resEMF;
  } else if (fJetType==PF) {
    delete [] fTChMult;     
    delete [] fTChHadFrac;  
    delete [] fTNeuHadFrac; 
    delete [] fTChEmFrac;   
    delete [] fTNeuEmFrac;  
    delete [] fTChMuFrac;  
    delete [] fTNeuMult;    	
    delete [] fTNConstituents;		
  }

}

//________________________________________________________________________________________
void JetFillerBase::createBranches(void) {

  addBranch("NJets",  "I", &fTnobj );
  addBranch("JPx",    "D", fTpx,   "NJets" );
  addBranch("JPy",    "D", fTpy,   "NJets" );
  addBranch("JPz",    "D", fTpz,   "NJets" );
  addBranch("JPt",    "D", fTpt,   "NJets" );
  addBranch("JE",     "D", fTe,    "NJets" );
  addBranch("JEt",    "D", fTet,   "NJets" );
  addBranch("JEta",   "D", fTeta,  "NJets" );
  addBranch("JPhi",   "D", fTphi,  "NJets" );
  addBranch("JScale", "D", fTscale,"NJets" );
  addBranch("JL1FastJetScale", "D", fTL1FastJetScale,"NJets" );
  addBranch("JArea",  "D", fTarea, "NJets" );
  addBranch("JbTagProbTkCntHighEff" , "D", fTjbTagProbTkCntHighEff , "NJets" );
  addBranch("JbTagProbTkCntHighPur" , "D", fTjbTagProbTkCntHighPur , "NJets" );
  addBranch("JbTagProbSimpSVHighEff", "D", fTjbTagProbSimpSVHighEff, "NJets" );
  addBranch("JbTagProbSimpSVHighPur", "D", fTjbTagProbSimpSVHighPur, "NJets" );
  addBranch("JIDLoose","I",fTIDLoose,"NJets" );

	
  if(fJetType==CALO) {
    addBranch("JID_HPD",        "D", fTID_HPD,       "NJets" );
    addBranch("JID_RBX",        "D", fTID_RBX,       "NJets" );
    addBranch("JID_n90Hits",    "D", fTID_n90Hits,   "NJets" );
    addBranch("JID_resEMF",     "D", fTID_resEMF,    "NJets" );   
    addBranch("JEMfrac",        "D", fTEMfrac,       "NJets" );
    addBranch("JNAssoTracks",   "I", fTjnAssoTracks, "NJets" );
    addBranch("JChfrac",        "D", fTjChfrac,      "NJets" );
    addBranch("JNConstituents", "I", fTNConstituents,"NJets" );				
    addBranch("Jn90"          , "I", fTn90,          "NJets" );				
  } else if (fJetType==JPT) {
    addBranch("JChMult",    "I", fTChMult,    "NJets");
    addBranch("JID_HPD",    "D",fTID_HPD     ,"NJets" );
    addBranch("JID_RBX",    "D",fTID_RBX     ,"NJets" );
    addBranch("JID_n90Hits","D",fTID_n90Hits ,"NJets" );
    addBranch("JID_resEMF", "D",fTID_resEMF  ,"NJets" );
  } else if (fJetType==PF) {
    addBranch("JChMult",    "I", fTChMult,    "NJets");   													
    addBranch("JNeuMult",   "I", fTNeuMult,    "NJets");		
    addBranch("JChHadfrac", "D", fTChHadFrac, "NJets" );
    addBranch("JNeuHadfrac","D", fTNeuHadFrac,"NJets" );
    addBranch("JChEmfrac",  "D", fTChEmFrac,  "NJets" );
    addBranch("JNeuEmfrac", "D", fTNeuEmFrac, "NJets" );
    addBranch("JChMufrac",  "D", fTChMuFrac,  "NJets" );
    addBranch("JNConstituents", "I", fTNConstituents,"NJets" );		
  }
}

//________________________________________________________________________________________
void JetFillerBase::reset(void) {

  fTnobj = 0;
	
  resetDouble(fTpx  ,gMaxnobjs);
  resetDouble(fTpy  ,gMaxnobjs);
  resetDouble(fTpz  ,gMaxnobjs);
  resetDouble(fTpt  ,gMaxnobjs);
  resetDouble(fTe   ,gMaxnobjs);
  resetDouble(fTet  ,gMaxnobjs);
  resetDouble(fTeta ,gMaxnobjs);
  resetDouble(fTphi ,gMaxnobjs);
  resetDouble(fTarea,gMaxnobjs);
  resetDouble(fTjbTagProbTkCntHighEff ,gMaxnobjs);
  resetDouble(fTjbTagProbTkCntHighPur ,gMaxnobjs);
  resetDouble(fTjbTagProbSimpSVHighEff,gMaxnobjs);
  resetDouble(fTjbTagProbSimpSVHighPur,gMaxnobjs);
  resetInt (fTIDLoose, gMaxnobjs);
	
  if (fJetType==CALO) {	
    resetDouble(fTID_RBX,gMaxnobjs);
    resetDouble(fTID_HPD,gMaxnobjs);
    resetDouble(fTID_n90Hits,gMaxnobjs);
    resetDouble(fTID_resEMF,gMaxnobjs);
    resetDouble(fTEMfrac,gMaxnobjs);
    resetInt(fTjnAssoTracks,gMaxnobjs);
    resetDouble(fTjChfrac, gMaxnobjs);
    resetInt(fTNConstituents,gMaxnobjs);		
  } else if (fJetType==JPT) {	
    resetDouble(fTID_RBX,gMaxnobjs);
    resetDouble(fTID_HPD,gMaxnobjs);
    resetDouble(fTID_n90Hits,gMaxnobjs);
    resetDouble(fTID_resEMF,gMaxnobjs);
    resetInt(fTChMult, gMaxnobjs);
  } else if (fJetType==PF) {
    resetInt(fTNeuMult, gMaxnobjs);
    resetInt(fTChMult, gMaxnobjs);
    resetDouble(fTChHadFrac,gMaxnobjs);
    resetDouble(fTNeuHadFrac,gMaxnobjs);
    resetDouble(fTChEmFrac,gMaxnobjs);
    resetDouble(fTNeuEmFrac,gMaxnobjs);
    resetDouble(fTChMuFrac,gMaxnobjs);
    resetInt(fTNConstituents,gMaxnobjs);
  }
}
