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
JetFillerReco::JetFillerReco( const edm::ParameterSet& config, TTree* tree, const bool& isRealData )
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


  edm::LogVerbatim("NTP") << " ==> JetFillerReco Constructor - " << fPrefix;
  edm::LogVerbatim("NTP") << "  Input Tag:        " << fTag.label();
  edm::LogVerbatim("NTP") << "---------------------------------";

}


//________________________________________________________________________________________
const int JetFillerReco::fillBranches(const edm::Event& iEvent, 
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
  Handle<JetTagCollection> jetsAndProbsTkCntHighEff;
  iEvent.getByLabel("trackCountingHighEffBJetTags",jetsAndProbsTkCntHighEff);

  Handle<JetTagCollection> jetsAndProbsTkCntHighPur;
  iEvent.getByLabel("trackCountingHighPurBJetTags",jetsAndProbsTkCntHighPur);

  Handle<JetTagCollection> jetsAndProbsSimpSVHighEff;
  iEvent.getByLabel("simpleSecondaryVertexHighEffBJetTags",jetsAndProbsSimpSVHighEff);

  Handle<JetTagCollection> jetsAndProbsSimpSVHighPur;
  iEvent.getByLabel("simpleSecondaryVertexHighPurBJetTags",jetsAndProbsSimpSVHighPur);

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
    fTpx[ijet]    = jet->px()*scale;
    fTpy[ijet]    = jet->py()*scale;
    fTpz[ijet]    = jet->pz()*scale;
    fTpt[ijet]    = jet->pt()*scale;
    fTe[ijet]     = jet->energy()*scale;
    fTet[ijet]    = jet->et()*scale;
    fTphi[ijet]   = jet->phi();
    fTeta[ijet]   = jet->eta();
    fTscale[ijet] = scale;

    // B-tagging probability (for 4 b-taggings)
    double mindr(999.99);
    for (unsigned int i = 0; i < jetsAndProbsTkCntHighEff->size(); i++){
      // Angular match between the two "collections"
      double deltar = reco::deltaR( jet->eta(), jet->phi(), 
                                    (*jetsAndProbsTkCntHighEff)[i].first->eta(), 
                                    (*jetsAndProbsTkCntHighEff)[i].first->phi());
      if( deltar <= fBtagMatchdeltaR && deltar < mindr)  {
        fTjbTagProbTkCntHighEff[ijet]=(*jetsAndProbsTkCntHighEff)[i].second;
        mindr = deltar;
      }
    }
    mindr = 999.99;
    for (unsigned int i = 0; i < jetsAndProbsTkCntHighPur->size(); i++){
      // Angular match between the two "collections"
      double deltar = reco::deltaR( jet->eta(), jet->phi(), 
                                    (*jetsAndProbsTkCntHighPur)[i].first->eta(), 
                                    (*jetsAndProbsTkCntHighPur)[i].first->phi());
      if( deltar <= fBtagMatchdeltaR && deltar < mindr)  {
        fTjbTagProbTkCntHighPur[ijet]=(*jetsAndProbsTkCntHighPur)[i].second;
        mindr = deltar;
      }
    }
    mindr = 999.99;
    for (unsigned int i = 0; i < jetsAndProbsSimpSVHighEff->size(); i++){
      // Angular match between the two "collections"
      double deltar = reco::deltaR( jet->eta(), jet->phi(), 
                                    (*jetsAndProbsSimpSVHighEff)[i].first->eta(), 
                                    (*jetsAndProbsSimpSVHighEff)[i].first->phi());
      if( deltar <= fBtagMatchdeltaR && deltar < mindr)  {
        fTjbTagProbSimpSVHighEff[ijet]=(*jetsAndProbsSimpSVHighEff)[i].second;
        mindr = deltar;
      }
    }
    mindr = 999.99;
    for (unsigned int i = 0; i < jetsAndProbsSimpSVHighPur->size(); i++){
      // Angular match between the two "collections"
      double deltar = reco::deltaR( jet->eta(), jet->phi(), 
                                    (*jetsAndProbsSimpSVHighPur)[i].first->eta(), 
                                    (*jetsAndProbsSimpSVHighPur)[i].first->phi());
      if( deltar <= fBtagMatchdeltaR && deltar < mindr)  {
        fTjbTagProbSimpSVHighPur[ijet]=(*jetsAndProbsSimpSVHighPur)[i].second;
        mindr = deltar;
      }
    }

    // -----------------------------------------
    // JPT jet specific
    if ( jetType()==JPT ) {
      reco::JetID jetID;
      const JPTJet* jptjet = static_cast<const JPTJet*>(&(*jet));
      edm::RefToBase<reco::Jet>  jetRef = jptjet->getCaloJetRef();
      jetID = (*jetsID)[ jetRef ];

      fTID_HPD[ijet]      = jetID.fHPD;
      fTID_RBX[ijet]      = jetID.fRBX;
      fTID_n90Hits[ijet]  = jetID.n90Hits;
      fTID_resEMF[ijet]   = jetID.restrictedEMF;

      fTChMult[ijet]      = jptjet->chargedMultiplicity();

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
      fTChHadFrac[ijet]     = CHF;
      fTNeuHadFrac[ijet]    = NHF;
      fTChEmFrac[ijet]      = CEF;
      fTNeuEmFrac[ijet]     = NEF;
      fTChMult[ijet]        = pjet->chargedMultiplicity();
      fTNeuMult[ijet]       = pjet->neutralMultiplicity(); 
      fTNConstituents[ijet] = pjet->nConstituents();  
    }

    // ------------------------------------------------------------
    // Calo jet specific
    if ( jetType()==CALO ) {
      reco::JetID jetID;
      const CaloJet* cjet = static_cast<const CaloJet*>(&(*jet));
      edm::RefToBase<reco::Jet> jetRef = jets->refAt(ijet);
      jetID = (*jetsID)[ jetRef ];

      fTID_HPD[ijet]      = jetID.fHPD;
      fTID_RBX[ijet]      = jetID.fRBX;
      fTID_n90Hits[ijet]  = jetID.n90Hits;
      fTID_resEMF[ijet]   = jetID.restrictedEMF;

      fTEMfrac[ijet]        = cjet->emEnergyFraction();
      fTNConstituents[ijet] = cjet->nConstituents();


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
      fTjnAssoTracks[ijet] = 0;
      fTjChfrac[ijet] = -1.; // Default (if jet-tracks association cone is outside tracker acceptance)
      if(fabs(jet->eta())<2.9){ // when the cone of dR=0.5 around the jet is (at least partially) inside the tracker acceptance
				// Tmp variables for vectorial sum of pt of tracks
        double pXtmp(0.), pYtmp(0.);
        // Loop over associated tracks:
        for(size_t t = 0; t < AssociatedTracks.size(); ++t){
          AssociatedTTracks.push_back(theB->build(AssociatedTracks[t])); // build transient tracks
          if(AssociatedTracks[t]->normalizedChi2()<10. && AssociatedTracks[t]->numberOfValidHits()>10 && AssociatedTracks[t]->pt()>1.){
            pXtmp += AssociatedTracks[t]->px();
            pYtmp += AssociatedTracks[t]->py();
            fTjnAssoTracks[ijet]++;
          }
        }
        fTjChfrac[ijet] = sqrt(pXtmp*pXtmp + pYtmp*pYtmp) / (jet->pt()*scale);

      } else { // The whole cone used for jet-tracks association is outside of the tracker acceptance
        fTjChfrac[ijet] = -1.;
      }
      AssociatedTracks.clear();
      AssociatedTTracks.clear();

    } // ----------------------------


    ++ijet;
  }
  fTnobj = ijet;

  return 0;

}
