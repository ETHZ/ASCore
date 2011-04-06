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
// $Id: LeptonFillerPat.h,v 1.3 2011/02/23 19:34:29 stiegerb Exp $
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

private:


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
  double* fTeta;
  double* fTphi;
  double* fTe;
  double* fTet;
  int*    fTcharge;

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
  fMinpt   = config.getParameter<double>("sel_minpt");
  fMaxeta  = config.getParameter<double>("sel_maxeta");
  gMaxnobjs = config.getUntrackedParameter<uint>("maxnobjs");
  
  fTag    = config.getUntrackedParameter<edm::InputTag>("tag");


  edm::LogVerbatim("NTP") << " ==> LeptonFillerPat Constructor - " << fPrefix;
  edm::LogVerbatim("NTP") << "  Input Tag:      " << fTag.label();
  edm::LogVerbatim("NTP") << "  Max n(objs):    " << gMaxnobjs;
  edm::LogVerbatim("NTP") << "  Min pt:         " << fMinpt;
  edm::LogVerbatim("NTP") << "  Max eta:        " << fMaxeta;
  edm::LogVerbatim("NTP") << "---------------------------------";

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
          
    fTparticleIso[pfqi]      = lepton.particleIso();
    fTchargedHadronIso[pfqi] = lepton.chargedHadronIso();
    fTneutralHadronIso[pfqi] = lepton.neutralHadronIso();
    fTphotonIso[pfqi]        = lepton.photonIso();
          
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

}

#endif
