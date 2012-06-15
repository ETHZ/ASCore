#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "DiLeptonAnalysis/NTupleProducer/interface/PFFiller.h"

#include <iostream>


//________________________________________________________________________________________
PFFiller::PFFiller( const edm::ParameterSet& config, const bool& isRealData )
  : FillerBase( config, isRealData )
{

  // Retrieve configuration parameters
  fPrefix          = config.getParameter<std::string>("prefix");
  fTag             = config.getParameter<edm::InputTag>("tag");
  fMinpt           = config.getParameter<double>("sel_minpt");
  fMaxeta          = config.getParameter<double>("sel_maxeta");

  gMaxnobjs = 200;

  edm::LogVerbatim("NTP") << " ==> PFFiller Constructor - " << fPrefix;
  edm::LogVerbatim("NTP") << "  Input Tag:        " << fTag.label();
  edm::LogVerbatim("NTP") << "---------------------------------";

}


//________________________________________________________________________________________
void PFFiller::fillProducts(edm::Event& iEvent,const edm::EventSetup& iSetup ) {

  using namespace edm;
  using namespace std;
  using namespace reco;

  // Retrieve collections
  Handle<View<reco::PFCandidate> > pfHandle;
  iEvent.getByLabel(fTag,pfHandle);
  View<reco::PFCandidate> candidates = *pfHandle;

  // PAT jets are already ordered by corrected pt: no need to re-order
  unsigned int icand(0);
  View<reco::PFCandidate>::const_iterator itend = candidates.end();
  for( View<reco::PFCandidate>::const_iterator it = candidates.begin();
       it != itend; ++it )
    {

      // Cut on corrected pT
      if(it->pt() < fMinpt) continue;

      // Save only the gMaxnjets first uncorrected jets
      if ( icand >= gMaxnobjs ) {
        edm::LogWarning("NTP") << "@SUB=FillBranches"
                               << "Maximum number of jets exceeded: "
                               << icand << " >= " << static_cast<int>(gMaxnobjs);
        break;
      }
      // Store the information (corrected)
      fTPx    ->push_back(it->px());
      fTPy    ->push_back(it->py());
      fTPz    ->push_back(it->pz());
      fTPt    ->push_back(it->pt());
      fTE     ->push_back(it->energy());
      fTEt    ->push_back(it->et());
      fTPhi   ->push_back(it->phi());
      fTEta   ->push_back(it->eta());
      fTVx    ->push_back(it->vx());
      fTVy    ->push_back(it->vy());
      fTVz    ->push_back(it->vz());

      // Types: https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideParticleFlow#ParticleTypes
      fTType  ->push_back(it->particleId());

      icand++;
    }
  *fTNObjs = icand;

}


//________________________________________________________________________________________
const std::vector<filler::PPair> PFFiller::declareProducts(void) {

  addProduct("NCandidates",   typeid(*fTNObjs));
  addProduct("Px",     typeid(*fTPx));
  addProduct("Py",     typeid(*fTPy));
  addProduct("Pz",     typeid(*fTPz));
  addProduct("Pt",     typeid(*fTPt));
  addProduct("E",      typeid(*fTE));
  addProduct("Et",     typeid(*fTEt));
  addProduct("Eta",    typeid(*fTEta));
  addProduct("Phi",    typeid(*fTPhi));
  addProduct("Type",   typeid(*fTType));
  addProduct("Vx",     typeid(*fTVx));
  addProduct("Vy",     typeid(*fTVy));
  addProduct("Vz",     typeid(*fTVz));

  return typeList;

}

//________________________________________________________________________________________
void PFFiller::resetProducts(void) {

  fTNObjs.reset(new int(0));
  
  fTPx     .reset(new std::vector<float>);
  fTPy     .reset(new std::vector<float>);
  fTPz     .reset(new std::vector<float>);
  fTPt     .reset(new std::vector<float>);
  fTE      .reset(new std::vector<float>);
  fTEt     .reset(new std::vector<float>);
  fTEta    .reset(new std::vector<float>);
  fTPhi    .reset(new std::vector<float>);
  fTType   .reset(new std::vector<int>);
  fTVx     .reset(new std::vector<float>);
  fTVy     .reset(new std::vector<float>);
  fTVz     .reset(new std::vector<float>);

}

//______________________________________________________________________________
void PFFiller::putProducts( edm::Event& e ) { 
  e.put(fTNObjs,fullName("NCandidates"));
	
  e.put(fTPx, fullName("Px"));
  e.put(fTPy, fullName("Py"));
  e.put(fTPz, fullName("Pz"));
  e.put(fTPt, fullName("Pt"));
  e.put(fTE,  fullName("E"));
  e.put(fTEt, fullName("Et"));
  e.put(fTEta,fullName("Eta"));
  e.put(fTPhi,fullName("Phi"));
  e.put(fTType,fullName("Type"));
  e.put(fTVx, fullName("Vx"));
  e.put(fTVy, fullName("Vy"));
  e.put(fTVz, fullName("Vz"));

}
