#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/JPTJetCollection.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/JetReco/interface/JetID.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/JetReco/interface/JetTracksAssociation.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"
#include "PhysicsTools/SelectorUtils/interface/JetIDSelectionFunctor.h"

#include "DiLeptonAnalysis/NTupleProducer/interface/JetFillerPat.h"

#include <iostream>


//________________________________________________________________________________________
JetFillerPat::JetFillerPat( const edm::ParameterSet& config, const bool& isRealData )
  : JetFillerBase( config, isRealData )
{

  // Retrieve configuration parameters
  fPrefix          = config.getParameter<std::string>("prefix");
  fTag             = config.getParameter<edm::InputTag>("tag");
  fMinpt           = config.getParameter<double>("sel_minpt");
  fMaxeta          = config.getParameter<double>("sel_maxeta");
  fBtagNames       = config.getParameter<std::vector<std::string> >("list_btags");

  // Check on number of allowed tags
  if ( fBtagNames.size()>gMaxNBtags )
    throw cms::Exception("BadConfig") << "Too many (" << fBtagNames.size() << ") btagging algos requested. "
                                      << "Maximum is " << gMaxNBtags;

  edm::LogVerbatim("NTP") << " ==> JetFillerPat Constructor - " << fPrefix;
  edm::LogVerbatim("NTP") << "  Input Tag:        " << fTag.label();
  edm::LogVerbatim("NTP") << "---------------------------------";

}


//________________________________________________________________________________________
void JetFillerPat::fillProducts(edm::Event& iEvent,const edm::EventSetup& iSetup ) {

  using namespace edm;
  using namespace std;
  using namespace reco;

  // Retrieve collections
  Handle<View<pat::Jet> > jetsHandle;
  iEvent.getByLabel(fTag,jetsHandle);
  View<pat::Jet> jets = *jetsHandle;

  // PFJetIDSelectionFunctor for LooseID.
  PFJetIDSelectionFunctor PFjetIDLoose( PFJetIDSelectionFunctor::FIRSTDATA, PFJetIDSelectionFunctor::LOOSE );
  JetIDSelectionFunctor   jetIDLoose( JetIDSelectionFunctor::PURE09, JetIDSelectionFunctor::LOOSE );

  // PAT jets are already ordered by corrected pt: no need to re-order
  unsigned int ijet(0);
  for( View<pat::Jet>::const_iterator Jit = jets.begin();
       Jit != jets.end(); ++Jit )
    {

      // Cut on corrected pT
      if(Jit->pt() < fMinpt) continue;

      // Save only the gMaxnjets first uncorrected jets
      if ( ijet >= gMaxnobjs ) {
        edm::LogWarning("NTP") << "@SUB=FillBranches"
                               << "Maximum number of jets exceeded: "
                               << ijet << " >= " << static_cast<int>(gMaxnobjs);
        break;
      }
      // Store the information (corrected)
      fTPx    ->push_back(Jit->px());
      fTPy    ->push_back(Jit->py());
      fTPz    ->push_back(Jit->pz());
      fTPt    ->push_back(Jit->pt());
      fTE     ->push_back(Jit->energy());
      fTEt    ->push_back(Jit->et());
      fTPhi   ->push_back(Jit->phi());
      fTEta   ->push_back(Jit->eta());
      fTArea  ->push_back(Jit->jetArea());

      fTScale ->push_back(1.0/Jit->jecFactor("Uncorrected")); // This is the inverse correction...
      fTL1FastJetScale ->push_back(Jit->jecFactor("L1FastJet")/Jit->jecFactor("Uncorrected")); 

      if(Jit->genJet())	fTFlavour ->push_back(Jit->partonFlavour());
      else fTFlavour ->push_back(0 );

      // B-tagging probability (for 8 b-taggings)
      size_t ibtag = 0;
      for ( std::vector<std::string>::const_iterator it = fBtagNames.begin();
         it != fBtagNames.end(); ++it ) {
         fTJbTagProb[ibtag]->push_back(Jit->bDiscriminator(*it));
         ++ibtag;
      }

      if(Jit->isPFJet()) fTIDLoose ->push_back((int) PFjetIDLoose(*Jit));
      else fTIDLoose ->push_back((int) jetIDLoose(*Jit));
      

      // -------------------------------------------------
      // PF jet specific
      if (fJetType==PF) {
       
	// energy fractions for JetID have to be computed w.r.t. uncorrected Jet!
	// all energy fractions (CHF, NHF, CEF, NEF) of pat::Jet apart from HFHadronEnergyFraction are anyway given w.r.t. uncorrected energy.
	// PFJet IDs calculated from these definitions of energy fractions are consistent with the result of the PFJetIDSelectionFunctor
	// moreover, the energy fractions add up to 1.
	pat::Jet const & jetuncorr=Jit->correctedJet("Uncorrected");

	// multiplicities do not depend on energy
        fTNConstituents ->push_back(Jit->numberOfDaughters()); // same as nConstituents
	fTChMult        ->push_back(Jit->chargedMultiplicity());
        fTNeuMult       ->push_back(Jit->neutralMultiplicity());

        fTChHadfrac  ->push_back(jetuncorr.chargedHadronEnergyFraction());
        fTNeuHadfrac ->push_back(jetuncorr.neutralHadronEnergyFraction() + jetuncorr.HFHadronEnergyFraction());
        fTChEmfrac   ->push_back(jetuncorr.chargedEmEnergyFraction());
        fTNeuEmfrac  ->push_back(jetuncorr.neutralEmEnergyFraction());
        fTChMufrac   ->push_back(jetuncorr.chargedMuEnergyFraction());    
	// Note: fTChHadfrac+fTNeuHadfrac+fTChEmfrac+fTNeuEmfrac+fTChMufrac==1

	fTPhofrac    ->push_back(jetuncorr.photonEnergyFraction()); // photons also count for neutralEmEnergy
        fTHFHadfrac  ->push_back(jetuncorr.HFHadronEnergyFraction());
        fTHFEMfrac   ->push_back(jetuncorr.HFEMEnergyFraction());   // also contained in neutralEmEnergy

      }

      ijet++;
    }
  *fTNObjs = ijet;

}
