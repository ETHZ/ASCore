#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/JPTJetCollection.h"
#include "DataFormats/JetReco/interface/JetID.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "DataFormats/JetReco/interface/JetTracksAssociation.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include "DiLeptonAnalysis/NTupleProducer/interface/JetFillerPat.h"
#include <iostream>


//________________________________________________________________________________________
JetFillerPat::JetFillerPat( const edm::ParameterSet& config, TTree* tree, 
                              const bool& isPat, const bool& isRealData )
  : JetFillerBase( config, tree, isPat, isRealData )
{
	
  // Retrieve configuration parameters
  fPrefix         = config.getUntrackedParameter<std::string>("prefix");
  fTag            = config.getUntrackedParameter<edm::InputTag>("tag");
  fMinpt          = config.getParameter<double>("sel_minpt");
  fMaxeta         = config.getParameter<double>("sel_maxeta");
  fJetCorrs       = config.getParameter<std::string>("corrections");
  fJetID          = config.getUntrackedParameter<edm::InputTag>("jet_id");		
  fJetTracksTag   = config.getUntrackedParameter<edm::InputTag>("tag_jetTracks");


  edm::LogVerbatim("NTP") << " ==> JetFillerPat Constructor - " << fPrefix;
  edm::LogVerbatim("NTP") << "  Input Tag:        " << fTag.label();
  edm::LogVerbatim("NTP") << "    fMinpt        = " << fMinpt;
  edm::LogVerbatim("NTP") << "    fMaxeta       = " << fMaxeta;
  edm::LogVerbatim("NTP") << "    fJetCorrs     = " << fJetCorrs;
  edm::LogVerbatim("NTP") << "    fJetID        = " << fJetID.label();
  edm::LogVerbatim("NTP") << "    fJetTracksTag = " << fJetTracksTag.label();
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
	  
//  const JetCorrector* jetCorr = JetCorrector::getJetCorrector(fJetCorrs,iSetup);
  
  // First loop: just get corrected pt and corresponding indices
//  unsigned int iraw(0);
//  vector<OrderPair> corrIndices; // Vector of indices and pt of corr. jets

  unsigned int ijet(0);

  for( View<pat::Jet>::const_iterator jet = jets.begin(); 
       jet != jets.end(); ++jet) //, ++iraw )
    {
    /*  // Store the (index,pt) pair, where pt is corrected
      Jet mijetraw = Jit->correctedJet(pat::JetCorrFactors::Raw);
      double scale = Jit->et()/mijetraw.et();
	//jetCorr->correction(Jit->p4());
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
    double scale = jetCorr->correction(jet->p4());
    if(jet->pt()*scale < fMinpt) continue;

    // Save only the gMaxnjets first uncorrected jets
    if (ijet >= gMaxnobjs){
      edm::LogWarning("NTP") << "@SUB=FillBranches"
                             << "Maximum number of jets exceeded: " 
                             << ijet << " >= " << static_cast<int>(gMaxnobjs);
      break;
    }
*/
    // Cut on corrected pT
    if(jet->pt() < fMinpt) continue;
    if (ijet >= gMaxnobjs){
      edm::LogWarning("NTP") << "@SUB=FillBranches"
                             << "Maximum number of jets exceeded: "
                             << ijet << " >= " << static_cast<int>(gMaxnobjs);
      break;
    }

    // Store the information (corrected)
    fTpx[ijet]    = jet->px();
    fTpy[ijet]    = jet->py();
    fTpz[ijet]    = jet->pz();
    fTpt[ijet]    = jet->pt();
    fTe[ijet]     = jet->energy();
    fTet[ijet]    = jet->et();
    fTphi[ijet]   = jet->phi();
    fTeta[ijet]   = jet->eta();

    Jet mijetraw = jet->correctedJet(pat::JetCorrFactors::Raw);
    double scale = jet->et()/mijetraw.et();
    fTscale[ijet] = scale;
		
    // -----------------------------------------
    // JPT jet specific
    if ( !(jet->isCaloJet()) && !(jet->isPFJet())){ // isJPTJet() ) { 38X.....
      reco::JetID jetID;
    //  const JPTJet* jptjet = static_cast<const JPTJet*>(&(*jet));
    //  edm::RefToBase<reco::Jet>  jetRef = jptjet->getCaloJetRef();

      //		jptjet->printJet();
			
      jetID =jet->jetID ();
			
      fTID_HPD[ijet]      = jetID.fHPD;
      fTID_RBX[ijet]      = jetID.fRBX;
      fTID_n90Hits[ijet]  = jetID.n90Hits;
      fTID_resEMF[ijet]   = jetID.restrictedEMF;
			
      fTChMult[ijet]      = jet->chargedMultiplicity(); 

    }	  
		
    // ------------------------------------------------------------
    // Calo jet specific
    if ( jet->isCaloJet()) {
      reco::JetID jetID;
//      const CaloJet* cjet = static_cast<const CaloJet*>(&(*jet));
 //     edm::RefToBase<reco::Jet> jetRef = jets->refAt(ijet);
      jetID = jet->jetID ();
				
      fTID_HPD[ijet]      = jetID.fHPD;
      fTID_RBX[ijet]      = jetID.fRBX;
      fTID_n90Hits[ijet]  = jetID.n90Hits;
      fTID_resEMF[ijet]   = jetID.restrictedEMF;
			
      fTEMfrac[ijet]        = jet->emEnergyFraction();
     std::vector<CaloTowerPtr> const & temp = jet->getCaloConstituents ();
      fTNConstituents[ijet] = temp.size(); //jet->n90();

			
      /////////////////////////////////////////////////////
      // calculate charge fraction
			
      // Jet-track association: get associated tracks
      vector<const reco::Track*> AssociatedTracks;
			
      const reco::TrackRefVector& tracks = jet->associatedTracks();
//JetTracksAssociation::getValue(*(jetTracksAssoc.product()),jetRef);
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
        fTjChfrac[ijet] = sqrt(pXtmp*pXtmp + pYtmp*pYtmp) / (jet->pt());
			
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
