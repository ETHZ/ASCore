#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"

#include "DiLeptonAnalysis/NTupleProducer/interface/JetFiller.h"


//________________________________________________________________________________________
JetFiller::JetFiller( const edm::ParameterSet& config, TTree* tree, 
                      const bool& isPat, const bool& isRealData )
  : fTree(tree),fIsPat(isPat),fIsRealData(isRealData)
{

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

  fTNConstituents = new int[gMaxnobjs];
  fTEMfrac     = new double[gMaxnobjs];
  fTHadFrac    = new double[gMaxnobjs];
  fTChHadFrac  = new double[gMaxnobjs];
  fTNeuHadFrac = new double[gMaxnobjs];
  fTChEmFrac   = new double[gMaxnobjs];
  fTNeuEmFrac  = new double[gMaxnobjs];
  fTMuonMultiplicity = new int[gMaxnobjs];
  
  
  // Retrieve configuration parameters
  fPrefix   = config.getUntrackedParameter<std::string>("prefix");
  fTag      = config.getUntrackedParameter<edm::InputTag>("tag");
  fMinpt    = config.getParameter<double>("sel_minpt");
  fMaxeta   = config.getParameter<double>("sel_maxeta");
  fJetCorrs = config.getParameter<std::string>("corrections");

  edm::LogVerbatim("NTP") << " ==> JetFiller Constructor - " << fPrefix;
  edm::LogVerbatim("NTP") << "  Input Tag: " << fTag.label();
  edm::LogVerbatim("NTP") << "    fMinpt      = " << fMinpt;
  edm::LogVerbatim("NTP") << "    fMaxeta     = " << fMaxeta;
  edm::LogVerbatim("NTP") << "    fJetCorrs   = " << fJetCorrs;
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
  delete [] fTNConstituents;
  delete [] fTEMfrac;
  delete [] fTHadFrac;
  delete [] fTChHadFrac;
  delete [] fTNeuHadFrac;
  delete [] fTChEmFrac;
  delete [] fTNeuEmFrac;
  delete [] fTMuonMultiplicity;
  
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
  addBranch("JNConstituents", "I", fTNConstituents,"NJets" );

  addBranch("JEMfrac", "D", fTEMfrac,"NJets" );
  addBranch("JHadfrac","D", fTHadFrac,"NJets" );

  addBranch("JChHadfrac", "D", fTChHadFrac, "NJets" );
  addBranch("JNeuHadfrac","D", fTNeuHadFrac,"NJets" );
  addBranch("JChEmfrac",  "D", fTChEmFrac,  "NJets" );
  addBranch("JNeuEmfrac", "D", fTNeuEmFrac, "NJets" );
  addBranch("JMuonMultiplicity","I", fTMuonMultiplicity,"NJets");

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
    fTNConstituents[ijet] = jet->nConstituents();
    
    //FIXME: NEED A NICE WAY TO DISTINGUISH JET SPECIES (PAT?)
    //      // Calo jet specific
    //      if ( jet->isCaloJet() ) {
    //         const CaloJet* jet = static_cast<const CaloJet*>(&(*jet));
    //         fTEMfrac[ijet]  = jet->emEnergyFraction();
    //         fTHadFrac[ijet] = jet->energyFractionHadronic();
    //      }
    //
    //      // PF jet specific
    //      if ( jet->isPFJet() ) {
    //         const PFJet* jet = static_cast<const PFJet*>(&(*jet));
    //         fTChHadFrac[ijet]  = jet->chargedHadronEnergyFraction();
    //         fTNeuHadFrac[ijet] = jet->neutralHadronEnergyFraction();
    //         fTChEmFrac[ijet]   = jet->chargedEmEnergyFraction();
    //         fTNeuEmFrac[ijet]  = jet->chargedEmEnergyFraction();
    //         fTMuonMultiplicity[ijet] = jet->muonMultiplicity();
    //      }
    
    ++ijet;
  }
  fTnobj = ijet;

  return 0;

}

//________________________________________________________________________________________
void JetFiller::reset(void) {

  fTnobj = 0;
  resetInt(fTMuonMultiplicity,gMaxnobjs);
  resetInt(fTNConstituents,gMaxnobjs);

  resetDouble(fTpx,gMaxnobjs);
  resetDouble(fTpy,gMaxnobjs);
  resetDouble(fTpz,gMaxnobjs);
  resetDouble(fTpt,gMaxnobjs);
  resetDouble(fTe,gMaxnobjs);
  resetDouble(fTet,gMaxnobjs);
  resetDouble(fTeta,gMaxnobjs);
  resetDouble(fTphi,gMaxnobjs);
  resetDouble(fTEMfrac,gMaxnobjs);
  resetDouble(fTHadFrac,gMaxnobjs);
  resetDouble(fTChHadFrac,gMaxnobjs);
  resetDouble(fTNeuHadFrac,gMaxnobjs);
  resetDouble(fTChEmFrac,gMaxnobjs);
  resetDouble(fTNeuEmFrac,gMaxnobjs);

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
