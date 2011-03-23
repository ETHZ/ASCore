#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/JPTJetCollection.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/JetReco/interface/JetID.h"
//#include "DataFormats/JetReco/interface/Jet.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/JetReco/interface/JetTracksAssociation.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

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
  fJetCorrs        = config.getParameter<std::string>("corrections");
  fJetID           = config.getUntrackedParameter<edm::InputTag>("jet_id");
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

  Handle<edm::ValueMap<reco::JetID> > jetsID;
  iEvent.getByLabel(fJetID,jetsID);

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
      fTscale[ijet] = 1.0/Jit->jecFactor("Uncorrected"); // This is the inverse correction...

      //      cout << "PFtoPAT jet " << Jit->p4() << " scale " << 1.0/Jit->jecFactor("Uncorrected") << endl;

      // B-tagging probability (for 4 b-taggings)
      fTjbTagProbTkCntHighEff [ijet] = Jit->bDiscriminator("trackCountingHighEffBJetTags"        );
      fTjbTagProbTkCntHighPur [ijet] = Jit->bDiscriminator("trackCountingHighPurBJetTags"        );
      fTjbTagProbSimpSVHighEff[ijet] = Jit->bDiscriminator("simpleSecondaryVertexHighEffBJetTags");
      fTjbTagProbSimpSVHighPur[ijet] = Jit->bDiscriminator("simpleSecondaryVertexHighPurBJetTags");



      // -------------------------------------------------
      // PF jet specific
      if (fJetType==PF) {
        double CHF=Jit->chargedHadronEnergyFraction();
        double NHF=Jit->neutralHadronEnergyFraction();
        double CEF=Jit->chargedEmEnergyFraction();
        double NEF=Jit->neutralEmEnergyFraction();
        double CMF=Jit->chargedMuEnergyFraction();
        
        double sum=CHF+NHF+CEF+NEF+CMF;
        if (sum >0) {
          CHF=CHF/sum;
          NHF=NHF/sum;
          CEF=CEF/sum;
          NEF=NEF/sum;
          CMF=CMF/sum;
        } else {
          edm::LogWarning("NTP") << "PFJets: energy fraction ==0 ";
        }
        fTChHadFrac[ijet]     = CHF;
        fTNeuHadFrac[ijet]    = NHF;
        fTChEmFrac[ijet]      = CEF;
        fTNeuEmFrac[ijet]     = NEF;
        fTChMult[ijet]        = Jit->chargedMultiplicity();
        fTNeuMult[ijet]       = Jit->neutralMultiplicity();
        fTNConstituents[ijet] = Jit->nConstituents();
      }

      ijet++;
    }
  fTnobj = ijet;

  return 0;

}
