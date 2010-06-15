#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/JPTJetCollection.h"
#include "DataFormats/JetReco/interface/JetID.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "DataFormats/JetReco/interface/JetTracksAssociation.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "DiLeptonAnalysis/NTupleProducer/interface/JetFiller.h"
#include <iostream>


//________________________________________________________________________________________
JetFiller::JetFiller( const edm::ParameterSet& config, TTree* tree, 
                      const bool& isPat, const bool& isRealData )
  : fTree(tree),fIsPat(isPat),fIsRealData(isRealData)
{
	
  // Retrieve configuration parameters
  fPrefix         = config.getUntrackedParameter<std::string>("prefix");
  fTag            = config.getUntrackedParameter<edm::InputTag>("tag");
  fJetID          = config.getUntrackedParameter<edm::InputTag>("jet_id");		
  fMinpt          = config.getParameter<double>("sel_minpt");
  fMaxeta         = config.getParameter<double>("sel_maxeta");
  fJetCorrs       = config.getParameter<std::string>("corrections");
  fJetTracksTag   = config.getUntrackedParameter<edm::InputTag>("tag_jetTracks");


  // parse fTag label for jet-type	
  std::string	jettype= fTag.label();
  if (std::string::npos != jettype.find("PF")){
	  fJetType = PF;
  }else if(std::string::npos != jettype.find("Calo")){
	  fJetType = CALO;  
  }else if(std::string::npos != jettype.find("JPT")){
	  fJetType = JPT;  	  
  }else{
		fJetType = unknown;
	  edm::LogWarning("NTP") << "!! Don't know JetType !!" << fTag.label();
  }	

  //FIXME: could be set from configuration file...
  gMaxnobjs = 100;

  // Define all arrays
  fTpx  = new double[gMaxnobjs];
  fTpy  = new double[gMaxnobjs];
  fTpz  = new double[gMaxnobjs];
  fTpt  = new double[gMaxnobjs];
  fTe   = new double[gMaxnobjs];
  fTet  = new double[gMaxnobjs];
  fTeta = new double[gMaxnobjs];
  fTphi = new double[gMaxnobjs];
  fTscale = new double[gMaxnobjs];
  
	if(fJetType==CALO){
		fTID_HPD       = new double[gMaxnobjs];
		fTID_RBX       = new double[gMaxnobjs];
		fTID_n90Hits   = new double[gMaxnobjs];
		fTID_resEMF    = new double[gMaxnobjs];
		fTEMfrac       = new double[gMaxnobjs];
		fTjChfrac      = new double[gMaxnobjs];
		fTjnAssoTracks = new int[gMaxnobjs];
		fTNConstituents = new int[gMaxnobjs];		
	}else if(fJetType==JPT){
		fTChMult     = new int[gMaxnobjs];
		fTID_HPD     = new double[gMaxnobjs];
	  fTID_RBX     = new double[gMaxnobjs];
	  fTID_n90Hits = new double[gMaxnobjs];
	  fTID_resEMF  = new double[gMaxnobjs];
	}else if(fJetType==PF){
		fTChMult     = new int[gMaxnobjs];
		fTChHadFrac  = new double[gMaxnobjs];
	  fTNeuHadFrac = new double[gMaxnobjs];
	  fTChEmFrac   = new double[gMaxnobjs];
	  fTNeuEmFrac  = new double[gMaxnobjs];
		fTNeuMult    = new int[gMaxnobjs];	
		fTNConstituents = new int[gMaxnobjs];
	}
	

  edm::LogVerbatim("NTP") << " ==> JetFiller Constructor - " << fPrefix;
  edm::LogVerbatim("NTP") << "  Input Tag:        " << fTag.label();
  edm::LogVerbatim("NTP") << "    fMinpt        = " << fMinpt;
  edm::LogVerbatim("NTP") << "    fMaxeta       = " << fMaxeta;
  edm::LogVerbatim("NTP") << "    fJetCorrs     = " << fJetCorrs;
  edm::LogVerbatim("NTP") << "    fJetID        = " << fJetID.label();
	edm::LogVerbatim("NTP") << "    fJetTracksTag = " << fJetTracksTag.label();
  edm::LogVerbatim("NTP") << "---------------------------------";

}

//________________________________________________________________________________________
JetFiller::~JetFiller(void) {
	
	
  // Delete all arrays
  delete [] fTpx;
  delete [] fTpy;
  delete [] fTpz;
  delete [] fTpt;
  delete [] fTe;
  delete [] fTet;
  delete [] fTeta;
  delete [] fTphi;
  delete [] fTscale;
	
	if(fJetType==CALO){
		delete [] fTID_HPD;     
		delete [] fTID_RBX;     
		delete [] fTID_n90Hits; 
		delete [] fTID_resEMF;  
		delete [] fTEMfrac;  
		delete [] fTjnAssoTracks;
		delete [] fTjChfrac;
		delete [] fTNConstituents;				
	}else if(fJetType==JPT){
		delete [] fTChMult;     
		delete [] fTID_HPD;    
		delete [] fTID_RBX;    
		delete [] fTID_n90Hits; 
		delete [] fTID_resEMF;
	}else if(fJetType==PF){
		delete [] fTChMult;     
		delete [] fTChHadFrac;  
	  delete [] fTNeuHadFrac; 
	  delete [] fTChEmFrac;   
	  delete [] fTNeuEmFrac;  
		delete [] fTNeuMult;    	
		delete [] fTNConstituents;		
	}

}

//________________________________________________________________________________________
void JetFiller::createBranches(void) {


  addBranch("NJets",  "I", &fTnobj );
  addBranch("JPx",    "D", fTpx, "NJets" );
  addBranch("JPy",    "D", fTpy, "NJets" );
  addBranch("JPz",    "D", fTpz, "NJets" );
  addBranch("JPt",    "D", fTpt, "NJets" );
  addBranch("JE",     "D", fTe,  "NJets" );
  addBranch("JEt",    "D", fTet, "NJets" );
  addBranch("JEta",   "D", fTeta,"NJets" );
  addBranch("JPhi",   "D", fTphi,"NJets" );
  addBranch("JScale", "D", fTscale,"NJets" );

	
	if(fJetType==CALO){
	  addBranch("JID_HPD",     "D",fTID_HPD,       "NJets" );
		addBranch("JID_RBX",     "D",fTID_RBX,       "NJets" );
	  addBranch("JID_n90Hits", "D",fTID_n90Hits,   "NJets" );
		addBranch("JID_resEMF",  "D",fTID_resEMF,    "NJets" );   
		addBranch("JEMfrac",     "D",fTEMfrac,       "NJets" );
		addBranch("JNAssoTracks","I",fTjnAssoTracks, "NJets" );
		addBranch("JChfrac",     "D",fTjChfrac,      "NJets" );
		addBranch("JNConstituents", "I", fTNConstituents,"NJets" );				
	}else if(fJetType==JPT){
		addBranch("JChMult",    "I", fTChMult,    "NJets");
	  addBranch("JID_HPD",    "D",fTID_HPD     ,"NJets" );
		addBranch("JID_RBX",    "D",fTID_RBX     ,"NJets" );
	  addBranch("JID_n90Hits","D",fTID_n90Hits ,"NJets" );
		addBranch("JID_resEMF", "D",fTID_resEMF  ,"NJets" );

	}else if(fJetType==PF){
		addBranch("JChMult",    "I", fTChMult,    "NJets");   													
		addBranch("JNeuMult",   "I", fTChMult,    "NJets");		
	  addBranch("JChHadfrac", "D", fTChHadFrac, "NJets" );
	  addBranch("JNeuHadfrac","D", fTNeuHadFrac,"NJets" );
	  addBranch("JChEmfrac",  "D", fTChEmFrac,  "NJets" );
	  addBranch("JNeuEmfrac", "D", fTNeuEmFrac, "NJets" );
		addBranch("JNConstituents", "I", fTNConstituents,"NJets" );		
	}
	
}

//________________________________________________________________________________________
const int JetFiller::fillBranches(const edm::Event& iEvent, 
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
	

  const JetCorrector* jetCorr = JetCorrector::getJetCorrector(fJetCorrs,iSetup);

  // First loop: just get corrected pt and corresponding indices
  unsigned int iraw(0);
  vector<OrderPair> corrIndices; // Vector of indices and pt of corr. jets
  for( View<Jet>::const_iterator Jit = jets->begin(); 
       Jit != jets->end(); ++Jit, ++iraw )
    {
      // Store the (index,pt) pair, where pt is corrected
      double scale = jetCorr->correction(Jit->p4());
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
		
		// -----------------------------------------
		// JPT jet specific
		if ( fJetType==JPT ) {
			reco::JetID jetID;
			const JPTJet* jptjet = static_cast<const JPTJet*>(&(*jet));
			edm::RefToBase<reco::Jet>  jetRef = jptjet->getCaloJetRef();
			
			//		jptjet->printJet();
			
			jetID = (*jetsID)[ jetRef ];
			
			fTID_HPD[ijet]      = jetID.fHPD;
			fTID_RBX[ijet]      = jetID.fRBX;
			fTID_n90Hits[ijet]  = jetID.n90Hits;
			fTID_resEMF[ijet]   = jetID.restrictedEMF;
			
			fTChMult[ijet]      = jptjet->chargedMultiplicity();

		}	  
		
		// -------------------------------------------------
		// PF jet specific
		if ( fJetType==PF ) {
			const PFJet* pjet = static_cast<const PFJet*>(&(*jet));
			
			double CHF=pjet->chargedHadronEnergyFraction();
			double NHF=pjet->neutralHadronEnergyFraction();
			double CEF=pjet->chargedEmEnergyFraction();
			double NEF=pjet->neutralEmEnergyFraction();
			double CMF=pjet->chargedMuEnergyFraction();
			
			double sum=CHF+NHF+CEF+NEF+CMF;
			if(sum >0){
				CHF=CHF/sum;
				NHF=NHF/sum;
				CEF=CEF/sum;
				NEF=NEF/sum;
				CMF=CMF/sum;
			}else{
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
    if ( fJetType==CALO ) {
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

//________________________________________________________________________________________
void JetFiller::reset(void) {

  fTnobj = 0;
	
  resetDouble(fTpx,gMaxnobjs);
  resetDouble(fTpy,gMaxnobjs);
  resetDouble(fTpz,gMaxnobjs);
  resetDouble(fTpt,gMaxnobjs);
  resetDouble(fTe,gMaxnobjs);
  resetDouble(fTet,gMaxnobjs);
  resetDouble(fTeta,gMaxnobjs);
  resetDouble(fTphi,gMaxnobjs);

	
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
		resetInt(fTNConstituents,gMaxnobjs);
  }	  

}


//________________________________________________________________________________________
void JetFiller::resetDouble(double* v, size_t size)
{
  for(size_t i = 0; i < size; ++i) v[i] = -999.9;
}

//________________________________________________________________________________________
void JetFiller::resetInt(int* v, size_t size)
{
  for(size_t i = 0; i < size; ++i) v[i] = -999;
}


//________________________________________________________________________________________
const bool JetFiller::addBranch(const char* name, const char* type, 
                                void* address, const char* size )
{
  

  // Form input
  std::string fullname(fPrefix+name);
  
  std::string branchType(fullname);
  if ( size ) // Size needs to be pre-fixed
    branchType += "[" + std::string(fPrefix+size) + "]";
  branchType += "/"+std::string(type);

  // Declare branch
  TBranch* b = fTree->Branch(fullname.c_str(),address,branchType.c_str());

  return !(b==0); // return 1 if branch was successfully created, 0 otherwise

}
