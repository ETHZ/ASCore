#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/JPTJetCollection.h"
#include "DataFormats/JetReco/interface/JetCollection.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/JetReco/interface/JetID.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/JetReco/interface/JetTracksAssociation.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include "DiLeptonAnalysis/NTupleProducer/interface/JetFillerReco.h"
#include <iostream>


//________________________________________________________________________________________
JetFillerReco::JetFillerReco( const edm::ParameterSet& config, const bool& isRealData )
  : JetFillerBase( config, isRealData )
{

  // Retrieve configuration parameters
  fPrefix          = config.getParameter<std::string>("prefix");
  fTag             = config.getParameter<edm::InputTag>("tag");
  fMinpt           = config.getParameter<double>("sel_minpt");
  fMaxeta          = config.getParameter<double>("sel_maxeta");
  fJetCorrs        = config.getParameter<std::string>("corrections");
  fJetID           = config.getParameter<edm::InputTag>("jet_id");		
  fJetTracksTag    = config.getParameter<edm::InputTag>("tag_jetTracks");


  edm::LogVerbatim("NTP") << " ==> JetFillerReco Constructor - " << fPrefix;
  edm::LogVerbatim("NTP") << "  Input Tag:        " << fTag.label();
  edm::LogVerbatim("NTP") << "---------------------------------";

}


//________________________________________________________________________________________
void JetFillerReco::fillProducts(edm::Event& iEvent, 
                                 const edm::EventSetup& iSetup ) {

  using namespace edm;
  using namespace std;
  using namespace reco;

  // Retrieve collections
  Handle<View<Jet> > jets;
  iEvent.getByLabel(fTag,jets);

  // Jet tracks association (already done in PAT)
  Handle<reco::JetTracksAssociation::Container> jetTracksAssoc;
  iEvent.getByLabel(fJetTracksTag,jetTracksAssoc);

  Handle<edm::ValueMap<reco::JetID> > jetsID;
  iEvent.getByLabel(fJetID,jetsID);

  // Get Transient Track Builder
  ESHandle<TransientTrackBuilder> theB;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);

  // collect information for b-tagging (4 tags)
  // FIXME: THIS SHOULD BE CONFIGURABLE!
  Handle<JetTagCollection> jetsAndProbsTkCntHighEff;
  Handle<JetTagCollection> jetsAndProbsTkCntHighPur;
  Handle<JetTagCollection> jetsAndProbsSimpSVHighEff;
  Handle<JetTagCollection> jetsAndProbsSimpSVHighPur;

  if( jetType()==PF ) {
    iEvent.getByLabel("newPFTrackCountingHighEffBJetTags",jetsAndProbsTkCntHighEff);
    iEvent.getByLabel("newPFTrackCountingHighPurBJetTags",jetsAndProbsTkCntHighPur);
    iEvent.getByLabel("newPFSimpleSecondaryVertexHighEffBJetTags",jetsAndProbsSimpSVHighEff);
    iEvent.getByLabel("newPFSimpleSecondaryVertexHighPurBJetTags",jetsAndProbsSimpSVHighPur);
  } else {
    iEvent.getByLabel("newTrackCountingHighEffBJetTags",jetsAndProbsTkCntHighEff);
    iEvent.getByLabel("newTrackCountingHighPurBJetTags",jetsAndProbsTkCntHighPur);
    iEvent.getByLabel("newSimpleSecondaryVertexHighEffBJetTags",jetsAndProbsSimpSVHighEff);
    iEvent.getByLabel("newSimpleSecondaryVertexHighPurBJetTags",jetsAndProbsSimpSVHighPur);
  }


  const JetCorrector* jetCorr = JetCorrector::getJetCorrector(fJetCorrs,iSetup);

  // First loop: just get corrected pt and corresponding indices
  unsigned int iraw(0);
  vector<OrderPair> corrIndices; // Vector of indices and pt of corr. jets
  for( View<Jet>::const_iterator Jit = jets->begin(); 
       Jit != jets->end(); ++Jit, ++iraw )
    {
      // Store the (index,pt) pair, where pt is corrected
      reco::JetBaseRef jetRef(edm::Ref<JetView>(jets,iraw));
      double scale = jetCorr->correction(*Jit,jetRef,iEvent,iSetup);
      corrIndices.push_back(make_pair(iraw,Jit->pt()*scale));
    }


  // Now sort indices by decreasing corrected pt
  IndexByPt indexComparator;
  std::sort(corrIndices.begin(), corrIndices.end(), indexComparator);

  // Second loop: get them ordered
  unsigned int ijet(0);
  for ( vector<OrderPair>::const_iterator it = corrIndices.begin(); 
	it != corrIndices.end(); ++it ) {

    unsigned int index = it->first;
    const Jet* jet = &((*jets)[index]);

    // Cut on corrected pT
    double scale = it->second/jet->pt();
    if(jet->pt()*scale < fMinpt) continue;

    // Save only the gMaxnjets first uncorrected jets
    if (ijet >= gMaxnobjs){
      edm::LogWarning("NTP") << "@SUB=FillBranches"
                             << "Maximum number of jets exceeded: " 
                             << ijet << " >= " << static_cast<int>(gMaxnobjs);
      break;
    }

    // Store the information (corrected)
    fTPx   ->push_back(jet->px()*scale);
    fTPy   ->push_back(jet->py()*scale);
    fTPz   ->push_back(jet->pz()*scale);
    fTPt   ->push_back(jet->pt()*scale);
    fTE    ->push_back(jet->energy()*scale);
    fTEt   ->push_back(jet->et()*scale);
    fTEta  ->push_back(jet->eta());
    fTPhi  ->push_back(jet->phi());
    fTScale->push_back(scale);

    // B-tagging probability (for 4 b-taggings)
    float btag_match_deltaR = 0.5;
    double mindr(999.99);


    // -----------------------------------------
    // JPT jet specific
    if ( jetType()==JPT ) {
      reco::JetID jetID;
      const JPTJet* jptjet = static_cast<const JPTJet*>(&(*jet));
      edm::RefToBase<reco::Jet>  jetRef = jptjet->getCaloJetRef();
      jetID = (*jetsID)[ jetRef ];

      fTID_HPD    ->push_back(jetID.fHPD);
      fTID_RBX    ->push_back(jetID.fRBX);
      fTID_n90Hits->push_back(jetID.n90Hits);
      fTID_resEMF ->push_back(jetID.restrictedEMF);
      fTChMult    ->push_back(jptjet->chargedMultiplicity());

    }	  

    // -------------------------------------------------
    // PF jet specific
    if (fJetType==PF) {
      const PFJet* pjet = static_cast<const PFJet*>(&(*jet));

      double CHF=pjet->chargedHadronEnergyFraction();
      double NHF=pjet->neutralHadronEnergyFraction();
      double CEF=pjet->chargedEmEnergyFraction();
      double NEF=pjet->neutralEmEnergyFraction();
      double CMF=pjet->chargedMuEnergyFraction();

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
      fTNConstituents->push_back(pjet->nConstituents());  
      fTChMult       ->push_back(pjet->chargedMultiplicity());
      fTNeuMult      ->push_back(pjet->neutralMultiplicity()); 
      fTChHadfrac    ->push_back(CHF);
      fTNeuHadfrac   ->push_back(NHF);
      fTChEmfrac     ->push_back(CEF);
      fTNeuEmfrac    ->push_back(NEF);
    }

    // ------------------------------------------------------------
    // Calo jet specific
    if ( jetType()==CALO ) {
      reco::JetID jetID;
      const CaloJet* cjet = static_cast<const CaloJet*>(&(*jet));
      edm::RefToBase<reco::Jet> jetRef = jets->refAt(ijet);
      jetID = (*jetsID)[ jetRef ];

      fTNConstituents->push_back(cjet->nConstituents());
      fTEMfrac       ->push_back(cjet->emEnergyFraction());
      fTID_HPD       ->push_back(jetID.fHPD);
      fTID_RBX       ->push_back(jetID.fRBX);
      fTID_n90Hits   ->push_back(jetID.n90Hits);
      fTn90          ->push_back(cjet->n90());
      fTID_resEMF    ->push_back(jetID.restrictedEMF);


      /////////////////////////////////////////////////////
      // calculate charge fraction

      // Jet-track association: get associated tracks
      vector<const reco::Track*> AssociatedTracks;

      const reco::TrackRefVector& tracks = JetTracksAssociation::getValue(*(jetTracksAssoc.product()),jetRef);
      for ( TrackRefVector::iterator it = tracks.begin(); it != tracks.end(); ++it ){
        AssociatedTracks.push_back( it->get() );
      }

      // Jet-track association: make transient tracks and store information
      vector<TransientTrack> AssociatedTTracks;
      fTNAssoTracks->push_back(0);
      fTChfrac->push_back(-1.); // Default (if jet-tracks association cone is outside tracker acceptance)
      if(fabs(jet->eta())<2.9) { // when the cone of dR=0.5 around the jet is (at least partially) inside the tracker acceptance
				// Tmp variables for vectorial sum of pt of tracks
        double pXtmp(0.), pYtmp(0.);
        // Loop over associated tracks:
        for(size_t t = 0; t < AssociatedTracks.size(); ++t){
          AssociatedTTracks.push_back(theB->build(AssociatedTracks[t])); // build transient tracks
          if(AssociatedTracks[t]->normalizedChi2()<10. && AssociatedTracks[t]->numberOfValidHits()>10 && AssociatedTracks[t]->pt()>1.){
            pXtmp += AssociatedTracks[t]->px();
            pYtmp += AssociatedTracks[t]->py();
            (*fTNAssoTracks)[ijet]++;
          }
        }
        (*fTChfrac)[ijet] = sqrt(pXtmp*pXtmp + pYtmp*pYtmp) / (jet->pt()*scale);
      }
      AssociatedTracks.clear();
      AssociatedTTracks.clear();

    } // ----------------------------


    ++ijet;
  }
  *fTNObjs = ijet;

}
