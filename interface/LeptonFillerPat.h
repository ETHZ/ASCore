#ifndef __DiLeptonAnalysis_NTupleProducer_LeptonFillerPat_H__
#define __DiLeptonAnalysis_NTupleProducer_LeptonFillerPat_H__
// 
// Package: NTupleProducer
// Class:   LeptonFillerPat
//
/* class LeptonFillerPat
   LeptonFillerPat.h
   Description:  generic class for basic jet dumper

*/
//
// $Id: LeptonFillerPat.h,v 1.3 2011/04/28 13:40:37 pnef Exp $
//
//

#include <string>
#include <vector>

#include "TTree.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"

#include "DiLeptonAnalysis/NTupleProducer/interface/FillerBase.h"

template <class LeptonType>
class LeptonFillerPat : public FillerBase {
public:
  /// Constructor: set pointer to tree
  LeptonFillerPat<LeptonType>( const edm::ParameterSet&, TTree* tree, 
                               const bool& isRealData );
  virtual ~LeptonFillerPat(void);


  /// Define all branches
  virtual void createBranches(void);
  /// Reset all branch containers
  virtual void reset(void);
  /// Fill all branches
  virtual const int fillBranches(const edm::Event&, const edm::EventSetup& );
  /// for specific leptons
  virtual void getSpecific(LeptonType lepton, size_t index);

  enum Type {
    El, Mu, Tau, unknown
  };
  Type fType;

private:


  /// Set and get jet type
  void setType( const Type& type ) { fType = type; }
  const Type getType(void) const { return fType; }
  
  //- Configuration parameters
  edm::InputTag fTag; 
  
  // Pre-selection
  double fMinpt;
  double fMaxeta;

  size_t gMaxnobjs;
  
  // Tree variables
  int     fTflagmaxexc;
  int     fTnobjstot;
  int     fTnobj;
  double* fTpx;
  double* fTpy;
  double* fTpz;
  double* fTpt;
  double* fTpterr;
  double* fTeta;
  double* fTphi;
  double* fTe;
  double* fTet;
  int*    fTcharge;
  int*    fTdecaymode;
  int*    fTelID80;
  int*    fTelID85;
  int*    fTelID90;
  int*    fTelID95;
  int*    fTmuNMatches;

  double* fTparticleIso;
  double* fTchargedHadronIso;
  double* fTneutralHadronIso;
  double* fTphotonIso;

};

typedef LeptonFillerPat<pat::Muon>     PatMuonFiller;
typedef LeptonFillerPat<pat::Electron> PatElectronFiller;
typedef LeptonFillerPat<pat::Tau>      PatTauFiller;


//________________________________________________________________________________________
template <class LeptonType>
LeptonFillerPat<LeptonType>::LeptonFillerPat( const edm::ParameterSet& config, TTree* tree, const bool& isRealData )
  : FillerBase( config, tree, isRealData )
{

  // Retrieve configuration parameters
  std::string leptontype    = config.getUntrackedParameter<std::string>("type");
  fMinpt                    = config.getParameter<double>("sel_minpt");
  fMaxeta                   = config.getParameter<double>("sel_maxeta");
  gMaxnobjs                 = config.getUntrackedParameter<uint>("maxnobjs");
  
  fTag                      = config.getUntrackedParameter<edm::InputTag>("tag");

  edm::LogVerbatim("NTP") << " ==> LeptonFillerPat Constructor - " << fPrefix;
  edm::LogVerbatim("NTP") << "  Input Tag:      " << fTag.label();
  edm::LogVerbatim("NTP") << "  Max n(objs):    " << gMaxnobjs;
  edm::LogVerbatim("NTP") << "  Min pt:         " << fMinpt;
  edm::LogVerbatim("NTP") << "  Max eta:        " << fMaxeta;
  edm::LogVerbatim("NTP") << "---------------------------------";

  if      (leptontype == "electron")   setType(El);
  else if (leptontype == "muon")       setType(Mu);
  else if (leptontype == "tau")        setType(Tau);  
  else {
    setType(unknown);
    edm::LogWarning("NTP") << "!! Don't know Lepton Type !!" << leptontype;
  }
  
  // Define all arrays
  fTpx                = new double[gMaxnobjs];
  fTpy                = new double[gMaxnobjs];
  fTpz                = new double[gMaxnobjs];
  fTpt                = new double[gMaxnobjs];
  fTeta               = new double[gMaxnobjs];
  fTphi               = new double[gMaxnobjs];
  fTe                 = new double[gMaxnobjs];
  fTet                = new double[gMaxnobjs];
  fTcharge            = new int[gMaxnobjs];

  fTparticleIso       = new double[gMaxnobjs];
  fTchargedHadronIso  = new double[gMaxnobjs];
  fTneutralHadronIso  = new double[gMaxnobjs];
  fTphotonIso         = new double[gMaxnobjs];

  if(fType == Tau){
  	fTdecaymode    = new int[gMaxnobjs];
  }else if(fType == El){
  	fTelID95       = new int[gMaxnobjs];
  	fTelID90       = new int[gMaxnobjs];
  	fTelID85       = new int[gMaxnobjs];
  	fTelID80       = new int[gMaxnobjs];
  }else if(fType == Mu){
  	fTmuNMatches   = new int[gMaxnobjs];
	fTpterr        = new double[gMaxnobjs];
  }

}


//________________________________________________________________________________________
template <class LeptonType>
const int LeptonFillerPat<LeptonType>::fillBranches(const edm::Event& iEvent,
                                        const edm::EventSetup& iSetup ) {

  // Retrieve collection
  edm::Handle<edm::View<LeptonType> > collection;
  iEvent.getByLabel(fTag,collection);

  size_t pfqi(0);  // Index of qualified leptons
  for (typename edm::View<LeptonType>::const_iterator it = collection->begin(); 
        it != collection->end(); ++it ) {
    // Check if maximum number of leptons is exceeded already:
    if(pfqi >= gMaxnobjs){
      edm::LogWarning("NTP") << "@SUB=analyze()"
                             << "Maximum number of " << fPrefix << " exceeded";
      fTflagmaxexc = 1;
      break;
    }
    fTnobjstot++;

    // PfMuon preselection:
    if (it->pt() < fMinpt) continue;
    if (fabs(it->eta()) > fMaxeta) continue;

    const LeptonType& lepton = (*it);
    fTpx[pfqi]     = lepton.px();
    fTpy[pfqi]     = lepton.py();
    fTpz[pfqi]     = lepton.pz();
    fTpt[pfqi]     = lepton.pt();
    fTeta[pfqi]    = lepton.eta();
    fTphi[pfqi]    = lepton.phi();
    fTe[pfqi]      = lepton.energy();
    fTet[pfqi]     = lepton.et();
    fTcharge[pfqi] = lepton.charge();
          
    fTparticleIso[pfqi]      = (lepton.chargedHadronIso()+lepton.neutralHadronIso()+lepton.photonIso())/lepton.pt();
    fTchargedHadronIso[pfqi] = lepton.chargedHadronIso();
    fTneutralHadronIso[pfqi] = lepton.neutralHadronIso();
    fTphotonIso[pfqi]        = lepton.photonIso();
  
    getSpecific(lepton, pfqi);
          
    ++pfqi;

  }
  fTnobj = pfqi;

  return 0;

}

//______________________________________________________________________________
template <class LeptonType>
LeptonFillerPat<LeptonType>::~LeptonFillerPat(void) {

  // Delete all arrays  
  delete [] fTpx;
  delete [] fTpy;
  delete [] fTpz;
  delete [] fTpt;
  delete [] fTe;
  delete [] fTet;
  delete [] fTeta;
  delete [] fTphi;
  delete [] fTcharge;

  delete [] fTparticleIso;
  delete [] fTchargedHadronIso;
  delete [] fTneutralHadronIso;
  delete [] fTphotonIso;

  if(fType == Tau){
  	delete [] fTdecaymode;
  }else if (fType == El){
  	delete [] fTelID95;
  	delete [] fTelID90;
  	delete [] fTelID85;
  	delete [] fTelID80;
  }else if (fType == Mu){
  	delete [] fTmuNMatches;
	delete [] fTpterr;
  }

}

//________________________________________________________________________________________
template <class LeptonType>
void LeptonFillerPat<LeptonType>::createBranches(void) {

  addBranch("MaxLepExc","I", &fTflagmaxexc);
  addBranch("NObjsTot", "I", &fTnobjstot );
  addBranch("NObjs",    "I", &fTnobj );
  addBranch("Px",     "D", fTpx,    "NObjs" );
  addBranch("Py",     "D", fTpy,    "NObjs" );
  addBranch("Pz",     "D", fTpz,    "NObjs" );
  addBranch("Pt",     "D", fTpt,    "NObjs" );
  addBranch("E",      "D", fTe,     "NObjs" );
  addBranch("Et",     "D", fTet,    "NObjs" );
  addBranch("Eta",    "D", fTeta,   "NObjs" );
  addBranch("Phi",    "D", fTphi,   "NObjs" );
  addBranch("Charge", "I", fTcharge,"NObjs" );

  addBranch("ParticleIso",      "D", fTparticleIso,  "NObjs" );
  addBranch("ChargedHadronIso", "D", fTchargedHadronIso,  "NObjs" );
  addBranch("NeutralHadronIso", "D", fTneutralHadronIso,  "NObjs" );
  addBranch("PhotonIso",        "D", fTphotonIso,  "NObjs" );
  
  if(fType == Tau){
  	addBranch("DecayMode", "I", fTdecaymode,"NObjs" );
  }else if(fType == El){
  	addBranch("ID95", "I", fTelID95, "NObjs");
  	addBranch("ID90", "I", fTelID90, "NObjs");
  	addBranch("ID85", "I", fTelID85, "NObjs");
  	addBranch("ID80", "I", fTelID80, "NObjs");
  }else if(fType == Mu){
  	addBranch("PtErr"   , "D", fTpterr      , "NObjs");
	addBranch("NMatches", "I", fTmuNMatches , "NObjs");
  }

}

//______________________________________________________________________________
template <class LeptonType>
void LeptonFillerPat<LeptonType>::reset(void) {

  fTflagmaxexc = 0;
  fTnobjstot = 0;
  fTnobj = 0;
  resetDouble(fTpx ,gMaxnobjs);
  resetDouble(fTpy ,gMaxnobjs);
  resetDouble(fTpz ,gMaxnobjs);
  resetDouble(fTpt ,gMaxnobjs);
  resetDouble(fTe  ,gMaxnobjs);
  resetDouble(fTet ,gMaxnobjs);
  resetDouble(fTeta,gMaxnobjs);
  resetDouble(fTphi,gMaxnobjs);
  resetInt(fTcharge,gMaxnobjs);

  resetDouble(fTparticleIso,gMaxnobjs);
  resetDouble(fTchargedHadronIso,gMaxnobjs);
  resetDouble(fTneutralHadronIso,gMaxnobjs);
  resetDouble(fTphotonIso,gMaxnobjs);

  if(fType==Tau){
  	resetInt(fTdecaymode,gMaxnobjs);
  }else if(fType ==El ){
  	resetInt(fTelID95, gMaxnobjs);
  	resetInt(fTelID90, gMaxnobjs);
  	resetInt(fTelID85, gMaxnobjs);
  	resetInt(fTelID80, gMaxnobjs);
  }else if(fType == Mu){
  	resetInt   (fTmuNMatches, gMaxnobjs);
	resetDouble(fTpterr     , gMaxnobjs);
  }
  

}

//________________________________________________________________________________________
template <class LeptonType>
void LeptonFillerPat<LeptonType>::getSpecific(LeptonType lepton, size_t index){
  	return;
}

//________________________________________________________________________________________
template <>
void LeptonFillerPat<pat::Tau>::getSpecific(pat::Tau lepton, size_t index){
	// speficic for PFTaus
	fTdecaymode[index]   = lepton.decayMode();	
	return;
}

template <>
void LeptonFillerPat<pat::Electron>::getSpecific(pat::Electron lepton, size_t index){
	// speficic for PFElectrons
	fTelID95[index]   = lepton.electronID("simpleEleId95cIso");
	fTelID90[index]   = lepton.electronID("simpleEleId90cIso");
	fTelID85[index]   = lepton.electronID("simpleEleId85cIso");
	fTelID80[index]   = lepton.electronID("simpleEleId80cIso");
	return;
}

template <>
void LeptonFillerPat<pat::Muon>::getSpecific(pat::Muon lepton, size_t index){
	// speficic for PFMuon
	fTmuNMatches[index]   = lepton.numberOfMatches();
	fTpterr     [index]   = lepton.globalTrack()->ptError();
	return;
}

#endif
