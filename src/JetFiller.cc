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
  gMaxnobjs = 20;

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

  edm::LogVerbatim("NTP") << "---------------------------------";
  edm::LogVerbatim("NTP") << " ==> JetFiller Constructor - " << fPrefix;
  edm::LogVerbatim("NTP") << "  Input Tag:" ;
  edm::LogVerbatim("NTP") << "    fTag        = " << fTag.label();

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

  unsigned int ijet(0);
  for( View<Jet>::const_iterator Jit = jets->begin(); 
       Jit != jets->end(); ++Jit )
    {
      // Cut on corrected pT
      const JetCorrector* jetCorr = JetCorrector::getJetCorrector(fJetCorrs,iSetup);
      double scale = jetCorr->correction(Jit->p4());
      if(Jit->pt()*scale < fMinpt) continue;

      // Save only the gMaxnjets first uncorrected jets
      if (ijet >= gMaxnobjs){
        edm::LogWarning("NTP") << "@SUB=FillBranches"
                               << "Maximum number of jets exceeded: " 
                               << ijet << " >= " << static_cast<int>(gMaxnobjs);
        break;
      }

      fTpx[ijet]  = Jit->px();
      fTpy[ijet]  = Jit->py();
      fTpz[ijet]  = Jit->pz();
      fTpt[ijet]  = Jit->pt();
      fTe[ijet]   = Jit->energy();
      fTet[ijet]  = Jit->et();
      fTphi[ijet] = Jit->phi();
      fTscale[ijet] = scale;
      fTNConstituents[ijet] = Jit->nConstituents();

      //FIXME: NEED A NICE WAY TO DISTINGUISH JET SPECIES (PAT?)
      //      // Calo jet specific
      //      if ( Jit->isCaloJet() ) {
      //         const CaloJet* jet = static_cast<const CaloJet*>(&(*Jit));
      //         fTEMfrac[ijet]  = jet->emEnergyFraction();
      //         fTHadFrac[ijet] = jet->energyFractionHadronic();
      //      }
      //
      //      // PF jet specific
      //      if ( Jit->isPFJet() ) {
      //         const PFJet* jet = static_cast<const PFJet*>(&(*Jit));
      //         fTChHadFrac[ijet]  = jet->chargedHadronEnergyFraction();
      //         fTNeuHadFrac[ijet] = jet->neutralHadronEnergyFraction();
      //         fTChEmFrac[ijet]   = jet->chargedEmEnergyFraction();
      //         fTNeuEmFrac[ijet]  = jet->chargedEmEnergyFraction();
      //         fTMuonMultiplicity[ijet] = jet->muonMultiplicity();
      //      }

      ++ijet;
    }
  fTnobj = ijet+1;

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
