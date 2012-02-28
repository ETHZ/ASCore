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
JetFillerPat::JetFillerPat( const edm::ParameterSet& config, TTree* tree, const bool& isRealData )
  : JetFillerBase( config, tree, isRealData )
{

  // Retrieve configuration parameters
  fPrefix          = config.getUntrackedParameter<std::string>("prefix");
  fTag             = config.getUntrackedParameter<edm::InputTag>("tag");
  fMinpt           = config.getParameter<double>("sel_minpt");
  fMaxeta          = config.getParameter<double>("sel_maxeta");
  fJetTracksTag    = config.getUntrackedParameter<edm::InputTag>("tag_jetTracks");
  fBtagMatchdeltaR = config.getParameter<double>("btag_matchdeltaR");


  edm::LogVerbatim("NTP") << " ==> JetFillerPat Constructor - " << fPrefix;
  edm::LogVerbatim("NTP") << "  Input Tag:        " << fTag.label();
  edm::LogVerbatim("NTP") << "---------------------------------";

}


//________________________________________________________________________________________
const int JetFillerPat::fillBranches(const edm::Event& iEvent,
                                      const edm::EventSetup& iSetup ) {

  using namespace edm;
  using namespace std;
  using namespace reco;

  // Retrieve collections
  Handle<View<pat::Jet> > jetsHandle;
  iEvent.getByLabel(fTag,jetsHandle);
  View<pat::Jet> jets = *jetsHandle;

  // Jet tracks association (already done in PAT)
  Handle<reco::JetTracksAssociation::Container> jetTracksAssoc;
  iEvent.getByLabel(fJetTracksTag,jetTracksAssoc);

  // Get Transient Track Builder
  ESHandle<TransientTrackBuilder> theB;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);

  // collect information for b-tagging (4 tags)
  Handle<JetTagCollection> jetsAndProbsTkCntHighEff;
  iEvent.getByLabel("trackCountingHighEffBJetTags",jetsAndProbsTkCntHighEff);

  Handle<JetTagCollection> jetsAndProbsTkCntHighPur;
  iEvent.getByLabel("trackCountingHighPurBJetTags",jetsAndProbsTkCntHighPur);

  Handle<JetTagCollection> jetsAndProbsSimpSVHighEff;
  iEvent.getByLabel("simpleSecondaryVertexHighEffBJetTags",jetsAndProbsSimpSVHighEff);

  Handle<JetTagCollection> jetsAndProbsSimpSVHighPur;
  iEvent.getByLabel("simpleSecondaryVertexHighPurBJetTags",jetsAndProbsSimpSVHighPur);

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
      fTpx[ijet]    = Jit->px();
      fTpy[ijet]    = Jit->py();
      fTpz[ijet]    = Jit->pz();
      fTpt[ijet]    = Jit->pt();
      fTe[ijet]     = Jit->energy();
      fTet[ijet]    = Jit->et();
      fTphi[ijet]   = Jit->phi();
      fTeta[ijet]   = Jit->eta();
      fTarea[ijet]  = Jit->jetArea();

      if(Jit->genJet()){
	fTflavour[ijet] = Jit->partonFlavour();
      }
      else
	fTflavour[ijet] = 0 ;

      fTscale[ijet] = 1.0/Jit->jecFactor("Uncorrected"); // This is the inverse correction...
      fTL1FastJetScale[ijet] = Jit->jecFactor("L1FastJet")/Jit->jecFactor("Uncorrected"); 

      // B-tagging probability (for 4 b-taggings)
      fTjbTagProbTkCntHighEff [ijet] = Jit->bDiscriminator("trackCountingHighEffBJetTags"        );
      fTjbTagProbTkCntHighPur [ijet] = Jit->bDiscriminator("trackCountingHighPurBJetTags"        );
      fTjbTagProbSimpSVHighEff[ijet] = Jit->bDiscriminator("simpleSecondaryVertexHighEffBJetTags");
      fTjbTagProbSimpSVHighPur[ijet] = Jit->bDiscriminator("simpleSecondaryVertexHighPurBJetTags");

      if(Jit->isPFJet()){
	fTIDLoose[ijet] = (int) PFjetIDLoose(*Jit);
      } else {
	fTIDLoose[ijet] = (int) jetIDLoose(*Jit);
      } 

      // -------------------------------------------------
      // PF jet specific
      if (fJetType==PF) {
       
	// energy fractions for JetID have to be computed w.r.t. uncorrected Jet!
	// all energy fractions (CHF, NHF, CEF, NEF) of pat::Jet apart from HFHadronEnergyFraction are anyway given w.r.t. uncorrected energy.
	// PFJet IDs calculated from these definitions of energy fractions are consistent with the result of the PFJetIDSelectionFunctor
	// moreover, the energy fractions add up to 1.
	pat::Jet const & jetuncorr=Jit->correctedJet("Uncorrected");
        fTChHadFrac[ijet]  =jetuncorr.chargedHadronEnergyFraction();
        fTNeuHadFrac[ijet] =jetuncorr.neutralHadronEnergyFraction() + jetuncorr.HFHadronEnergyFraction();
        fTChEmFrac[ijet]   =jetuncorr.chargedEmEnergyFraction();
        fTNeuEmFrac[ijet]  =jetuncorr.neutralEmEnergyFraction();
        fTChMuFrac[ijet]   =jetuncorr.chargedMuEnergyFraction();
	fTPhoFrac  [ijet] = jetuncorr.photonEnergyFraction(); // photons also count for neutralEmEnergy
	fTHFHadFrac[ijet] = jetuncorr.HFHadronEnergyFraction();
	fTHFEMFrac [ijet] = jetuncorr.HFEMEnergyFraction();   // also contained in neutralEmEnergy
    
	// multiplicities do not depend on energy
	fTChMult[ijet]        = Jit->chargedMultiplicity();
        fTNeuMult[ijet]       = Jit->neutralMultiplicity();
        fTNConstituents[ijet] = Jit->numberOfDaughters(); // same as nConstituents
    
      }

      ijet++;
    }
  fTnobj = ijet;

  return 0;

}
