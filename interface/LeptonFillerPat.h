#ifndef __DiLeptonAnalysis_NTupleProducer_LeptonFillerPat_H__
#define __DiLeptonAnalysis_NTupleProducer_LeptonFillerPat_H__
// 
// Package: NTupleProducer
// Class:   LeptonFillerPat
//
/* class LeptonFillerPat
   LeptonFillerPat.h
   Description:  generic class for basic jet dumper

*/
//
// $Id: LeptonFillerPat.h,v 1.5.2.2 2012/01/27 15:07:22 fronga Exp $
//
//

#include <string>
#include <vector>

#include "TTree.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"

#include "DiLeptonAnalysis/NTupleProducer/interface/FillerBase.h"

template <class LeptonType>
class LeptonFillerPat : public FillerBase {
public:
    /// Constructor: set pointer to tree
    LeptonFillerPat<LeptonType>( const edm::ParameterSet&, const bool& isRealData );
    virtual ~LeptonFillerPat(void) {}

    /// Define all branches
    virtual const std::vector<filler::PPair> declareProducts(void);
    /// Reset all branch containers
    virtual void resetProducts(void);
    /// Fill all branches
    virtual void fillProducts(edm::Event&, const edm::EventSetup& );
    /// Put products in the event data
    virtual void putProducts( edm::Event& );

    enum Type {
        El, Mu, Tau, unknown
    };
    Type fType;

private:


    /// retrieve specific lepton information
    virtual void getSpecific(const LeptonType& lepton) {}

    /// Set and get jet type
    void setType( const Type& type ) { fType = type; }
    const Type getType(void) const { return fType; }
  
    //- Configuration parameters
    edm::InputTag fTag; 
  
    // Pre-selection
    double fMinpt;
    double fMaxeta;

    size_t gMaxnobjs;
  
    // Tree variables
    std::auto_ptr<int>     fTMaxLepExc;
    std::auto_ptr<int>     fTNObjsTot;
    std::auto_ptr<int>     fTNObjs;
    std::auto_ptr<std::vector<float> >  fTPx;
    std::auto_ptr<std::vector<float> >  fTPy;
    std::auto_ptr<std::vector<float> >  fTPz;
    std::auto_ptr<std::vector<float> >  fTPt;
    std::auto_ptr<std::vector<float> >  fTPtErr;
    std::auto_ptr<std::vector<float> >  fTEta;
    std::auto_ptr<std::vector<float> >  fTPhi;
    std::auto_ptr<std::vector<float> >  fTE;
    std::auto_ptr<std::vector<float> >  fTEt;
    std::auto_ptr<std::vector<int> >    fTCharge;
    std::auto_ptr<std::vector<int> >    fTDecayMode;
    std::auto_ptr<std::vector<float> >  fTVz; 
    std::auto_ptr<std::vector<float> >  fTEmFraction; 
    std::auto_ptr<std::vector<float> >  fTJetPt;
    std::auto_ptr<std::vector<float> >  fTJetEta;
    std::auto_ptr<std::vector<float> >  fTJetPhi;
    std::auto_ptr<std::vector<float> >  fTJetMass;
    std::auto_ptr<std::vector<float> >  fTLeadingTkPt;
    std::auto_ptr<std::vector<float> >  fTLeadingNeuPt;
    std::auto_ptr<std::vector<float> >  fTLeadingTkHcalenergy;
    std::auto_ptr<std::vector<float> >  fTLeadingTkEcalenergy;
    std::auto_ptr<std::vector<int> >    fTNumChargedHadronsSignalCone;
    std::auto_ptr<std::vector<int> >    fTNumNeutralHadronsSignalCone;
    std::auto_ptr<std::vector<int> >    fTNumPhotonsSignalCone;
    std::auto_ptr<std::vector<int> >    fTNumParticlesSignalCone;
    std::auto_ptr<std::vector<int> >    fTNumChargedHadronsIsoCone;
    std::auto_ptr<std::vector<int> >    fTNumNeutralHadronsIsoCone;
    std::auto_ptr<std::vector<int> >    fTNumPhotonsIsolationCone;
    std::auto_ptr<std::vector<int> >    fTNumParticlesIsolationCone;
    std::auto_ptr<std::vector<float> >  fTPtSumChargedParticlesIsoCone;
    std::auto_ptr<std::vector<float> >  fTPtSumPhotonsIsoCone;
    std::auto_ptr<std::vector<float> >  fTDecayModeFinding;
    std::auto_ptr<std::vector<float> >  fTVLooseIso;
    std::auto_ptr<std::vector<float> >  fTLooseIso;
    std::auto_ptr<std::vector<float> >  fTTightIso;
    std::auto_ptr<std::vector<float> >  fTMediumIso;
/*   std::auto_ptr<std::vector<float> >  fTVLooseChargedIso; */
/*   std::auto_ptr<std::vector<float> >  fTLooseChargedIso; */
/*   std::auto_ptr<std::vector<float> >  fTTightChargedIso; */
/*   std::auto_ptr<std::vector<float> >  fTMediumChargedIso; */
/*   std::auto_ptr<std::vector<float> >  fTVLooseIsoDBSumPtCorr; */
/*   std::auto_ptr<std::vector<float> >  fTLooseIsoDBSumPtCorr; */
/*   std::auto_ptr<std::vector<float> >  fTTightIsoDBSumPtCorr; */
/*   std::auto_ptr<std::vector<float> >  fTMediumIsoDBSumPtCorr; */
/*   std::auto_ptr<std::vector<float> >  fTVLooseCombinedIsoDBSumPtCorr; */
/*   std::auto_ptr<std::vector<float> >  fTLooseCombinedIsoDBSumPtCorr; */
/*   std::auto_ptr<std::vector<float> >  fTTightCombinedIsoDBSumPtCorr; */
/*   std::auto_ptr<std::vector<float> >  fTMediumCombinedIsoDBSumPtCorr; */
    std::auto_ptr<std::vector<float> >  fTLooseElectronRejection;
    std::auto_ptr<std::vector<float> >  fTTightElectronRejection;
    std::auto_ptr<std::vector<float> >  fTMediumElectronRejection;
    std::auto_ptr<std::vector<float> >  fTLooseMuonRejection;
    std::auto_ptr<std::vector<float> >  fTTightMuonRejection;


    std::auto_ptr<std::vector<int> >    fTID80;
    std::auto_ptr<std::vector<int> >    fTID85;
    std::auto_ptr<std::vector<int> >    fTID90;
    std::auto_ptr<std::vector<int> >    fTID95;
    std::auto_ptr<std::vector<int> >    fTMuNMatches;

    std::auto_ptr<std::vector<float> > fTParticleIso;
    std::auto_ptr<std::vector<float> > fTChargedHadronIso;
    std::auto_ptr<std::vector<float> > fTNeutralHadronIso;
    std::auto_ptr<std::vector<float> > fTPhotonIso;

};

typedef LeptonFillerPat<pat::Muon>     PatMuonFiller;
typedef LeptonFillerPat<pat::Electron> PatElectronFiller;
typedef LeptonFillerPat<pat::Tau>      PatTauFiller;


//________________________________________________________________________________________
template <class LeptonType>
LeptonFillerPat<LeptonType>::LeptonFillerPat( const edm::ParameterSet& config, const bool& isRealData )
    : FillerBase( config, isRealData )
{

    // Retrieve configuration parameters
    std::string leptontype    = config.getParameter<std::string>("type");
    fMinpt                    = config.getParameter<double>("sel_minpt");
    fMaxeta                   = config.getParameter<double>("sel_maxeta");
    gMaxnobjs                 = config.getParameter<uint>("maxnobjs");
  
    fTag                      = config.getParameter<edm::InputTag>("tag");

    edm::LogVerbatim("NTP") << " ==> LeptonFillerPat Constructor - " << fPrefix;
    edm::LogVerbatim("NTP") << "  Input Tag:      " << fTag.label();
    edm::LogVerbatim("NTP") << "  Max n(objs):    " << gMaxnobjs;
    edm::LogVerbatim("NTP") << "  Min pt:         " << fMinpt;
    edm::LogVerbatim("NTP") << "  Max eta:        " << fMaxeta;
    edm::LogVerbatim("NTP") << "---------------------------------";

    if      (leptontype == "electron")   setType(El);
    else if (leptontype == "muon")       setType(Mu);
    else if (leptontype == "tau")        setType(Tau);  
    else {
        setType(unknown);
        edm::LogWarning("NTP") << "!! Don't know Lepton Type !!" << leptontype;
    }
  
}


//________________________________________________________________________________________
template <class LeptonType>
void LeptonFillerPat<LeptonType>::fillProducts(edm::Event& iEvent,
                                              const edm::EventSetup& iSetup ) {

    // Retrieve collection
    edm::Handle<edm::View<LeptonType> > collection;
    iEvent.getByLabel(fTag,collection);
    
    size_t pfqi(0);  // Index of qualified leptons
    for (typename edm::View<LeptonType>::const_iterator it = collection->begin(); 
         it != collection->end(); ++it ) {
        // Check if maximum number of leptons is exceeded already:
        if(pfqi >= gMaxnobjs){
            edm::LogWarning("NTP") << "@SUB=analyze()"
                                   << "Maximum number of " << fPrefix << " leptons exceeded";
            *fTMaxLepExc = 1;
            break;
        }
        (*fTNObjsTot)++;

        // PfMuon preselection:
        if (it->pt() < fMinpt) continue;
        if (fabs(it->eta()) > fMaxeta) continue;
        
        const LeptonType& lepton = (*it);
        fTPx    ->push_back( lepton.px() );
        fTPy    ->push_back( lepton.py() );
        fTPz    ->push_back( lepton.pz() );
        fTPt    ->push_back( lepton.pt() );
        fTEta   ->push_back( lepton.eta() );
        fTPhi   ->push_back( lepton.phi() );
        fTE     ->push_back( lepton.energy() );
        fTEt    ->push_back( lepton.et() );
        fTCharge->push_back( lepton.charge() );
          
        fTParticleIso     ->push_back( (lepton.chargedHadronIso()+lepton.neutralHadronIso()+lepton.photonIso())/lepton.pt() );
        fTChargedHadronIso->push_back( lepton.chargedHadronIso() );
        fTNeutralHadronIso->push_back( lepton.neutralHadronIso() );
        fTPhotonIso       ->push_back( lepton.photonIso() );
  
        getSpecific(lepton);
          
        ++pfqi;

    }
    *fTNObjs = pfqi;

}

//________________________________________________________________________________________
template <class LeptonType>
const std::vector<filler::PPair> LeptonFillerPat<LeptonType>::declareProducts(void) {

    addProduct("MaxLepExc",typeid(*fTMaxLepExc));
    addProduct("NObjsTot", typeid(*fTNObjsTot));
    addProduct("NObjs",    typeid(*fTNObjs));

    addProduct("Px",       typeid(*fTPx));
    addProduct("Py",       typeid(*fTPy));
    addProduct("Pz",       typeid(*fTPz));
    addProduct("Pt",       typeid(*fTPt));
    addProduct("E",        typeid(*fTE));
    addProduct("Et",       typeid(*fTEt));
    addProduct("Eta",      typeid(*fTEta));
    addProduct("Phi",      typeid(*fTPhi));
    addProduct("Charge",   typeid(*fTCharge));

    addProduct("ParticleIso",      typeid(*fTParticleIso));
    addProduct("ChargedHadronIso", typeid(*fTChargedHadronIso));
    addProduct("NeutralHadronIso", typeid(*fTNeutralHadronIso));
    addProduct("PhotonIso",        typeid(*fTPhotonIso));
  
    if(fType == Tau){
  	addProduct("DecayMode", typeid(*fTDecayMode));
  	addProduct("Vz",        typeid(*fTVz)); 
  	addProduct("EmFraction",typeid(*fTEmFraction)); 
  	addProduct("JetPt",     typeid(*fTJetPt));
  	addProduct("JetEta",    typeid(*fTJetEta));
  	addProduct("JetPhi",    typeid(*fTJetPhi));
  	addProduct("JetMass",   typeid(*fTJetMass));
  	addProduct("LeadingTkPt", typeid(*fTLeadingTkPt));
  	addProduct("LeadingNeuPt",typeid(*fTLeadingNeuPt));
  	addProduct("LeadingTkHcalenergy", typeid(*fTLeadingTkHcalenergy));
  	addProduct("LeadingTkEcalenergy", typeid(*fTLeadingTkEcalenergy));
  	addProduct("NumChargedHadronsSignalCone", typeid(*fTNumChargedHadronsSignalCone));
  	addProduct("NumNeutralHadronsSignalCone", typeid(*fTNumNeutralHadronsSignalCone));
  	addProduct("NumPhotonsSignalCone",     typeid(*fTNumPhotonsSignalCone));
  	addProduct("NumParticlesSignalCone",   typeid(*fTNumParticlesSignalCone));
  	addProduct("NumChargedHadronsIsoCone", typeid(*fTNumChargedHadronsIsoCone));
  	addProduct("NumNeutralHadronsIsoCone", typeid(*fTNumNeutralHadronsIsoCone));
  	addProduct("NumPhotonsIsolationCone",  typeid(*fTNumPhotonsIsolationCone));
  	addProduct("NumParticlesIsolationCone",typeid(*fTNumParticlesIsolationCone));
  	addProduct("PtSumChargedParticlesIsoCone", typeid(*fTPtSumChargedParticlesIsoCone));
  	addProduct("PtSumPhotonsIsoCone",      typeid(*fTPtSumPhotonsIsoCone));
	addProduct("DecayModeFinding", typeid(*fTDecayModeFinding));
  	addProduct("VLooseIso", typeid(*fTVLooseIso));
  	addProduct("LooseIso",  typeid(*fTLooseIso));
  	addProduct("TightIso",  typeid(*fTTightIso));
  	addProduct("MediumIso", typeid(*fTMediumIso));
/*   	addProduct("VLooseChargedIso", typeid(*fTVLooseChargedIso)); */
/*   	addProduct("LooseChargedIso",  typeid(*fTLooseChargedIso)); */
/*   	addProduct("TightChargedIso",  typeid(*fTTightChargedIso)); */
/*   	addProduct("MediumChargedIso", typeid(*fTMediumChargedIso)); */
/*   	addProduct("VLooseIsoDBSumPtCorr", typeid(*fTVLooseIsoDBSumPtCorr)); */
/*   	addProduct("LooseIsoDBSumPtCorr",  typeid(*fTLooseIsoDBSumPtCorr)); */
/*   	addProduct("TightIsoDBSumPtCorr",  typeid(*fTTightIsoDBSumPtCorr)); */
/*   	addProduct("MediumIsoDBSumPtCorr", typeid(*fTMediumIsoDBSumPtCorr)); */
/*   	addProduct("VLooseCombinedIsoDBSumPtCorr", typeid(*fTVLooseCombinedIsoDBSumPtCorr)); */
/*   	addProduct("LooseCombinedIsoDBSumPtCorr",  typeid(*fTLooseCombinedIsoDBSumPtCorr)); */
/*   	addProduct("TightCombinedIsoDBSumPtCorr",  typeid(*fTTightCombinedIsoDBSumPtCorr)); */
/*   	addProduct("MediumCombinedIsoDBSumPtCorr", typeid(*fTMediumCombinedIsoDBSumPtCorr)); */
  	addProduct("LooseElectronRejection", typeid(*fTLooseElectronRejection));
  	addProduct("TightElectronRejection", typeid(*fTTightElectronRejection));
  	addProduct("MediumElectronRejection",typeid(*fTMediumElectronRejection));
  	addProduct("LooseMuonRejection", typeid(*fTLooseMuonRejection));
  	addProduct("TightMuonRejection", typeid(*fTTightMuonRejection));
    }else if(fType == El){
  	addProduct("ID95", typeid(*fTID95));
  	addProduct("ID90", typeid(*fTID90));
  	addProduct("ID85", typeid(*fTID85));
  	addProduct("ID80", typeid(*fTID80));
    }else if(fType == Mu){
  	addProduct("PtErr"   , typeid(*fTPtErr));
	addProduct("NMatches", typeid(*fTMuNMatches));
    }

    return typeList;

}

//________________________________________________________________________________________
template <class LeptonType>
void LeptonFillerPat<LeptonType>::putProducts(edm::Event& e) {

    e.put(fTMaxLepExc,fullName("MaxLepExc"));
    e.put(fTNObjsTot,fullName("NObjsTot"));
    e.put(fTNObjs,fullName("NObjs"));

    e.put(fTPx,fullName("Px"));
    e.put(fTPy,fullName("Py"));
    e.put(fTPz,fullName("Pz"));
    e.put(fTPt,fullName("Pt"));
    e.put(fTE,fullName("E"));
    e.put(fTEt,fullName("Et"));
    e.put(fTEta,fullName("Eta"));
    e.put(fTPhi,fullName("Phi"));
    e.put(fTCharge,fullName("Charge"));

    e.put(fTParticleIso,fullName("ParticleIso"));
    e.put(fTChargedHadronIso,fullName("ChargedHadronIso"));
    e.put(fTNeutralHadronIso,fullName("NeutralHadronIso"));
    e.put(fTPhotonIso,fullName("PhotonIso"));
  
    if(fType == Tau){
  	e.put(fTDecayMode,fullName("DecayMode"));
  	e.put(fTVz,fullName("Vz")); 
  	e.put(fTEmFraction,fullName("EmFraction")); 
  	e.put(fTJetPt,fullName("JetPt"));
  	e.put(fTJetEta,fullName("JetEta"));
  	e.put(fTJetPhi,fullName("JetPhi"));
  	e.put(fTJetMass,fullName("JetMass"));
  	e.put(fTLeadingTkPt,fullName("LeadingTkPt"));
  	e.put(fTLeadingNeuPt,fullName("LeadingNeuPt"));
  	e.put(fTLeadingTkHcalenergy,fullName("LeadingTkHcalenergy"));
  	e.put(fTLeadingTkEcalenergy,fullName("LeadingTkEcalenergy"));
  	e.put(fTNumChargedHadronsSignalCone,fullName("NumChargedHadronsSignalCone"));
  	e.put(fTNumNeutralHadronsSignalCone,fullName("NumNeutralHadronsSignalCone"));
  	e.put(fTNumPhotonsSignalCone,fullName("NumPhotonsSignalCone"));
  	e.put(fTNumParticlesSignalCone,fullName("NumParticlesSignalCone"));
  	e.put(fTNumChargedHadronsIsoCone,fullName("NumChargedHadronsIsoCone"));
  	e.put(fTNumNeutralHadronsIsoCone,fullName("NumNeutralHadronsIsoCone"));
  	e.put(fTNumPhotonsIsolationCone,fullName("NumPhotonsIsolationCone"));
  	e.put(fTNumParticlesIsolationCone,fullName("NumParticlesIsolationCone"));
  	e.put(fTPtSumChargedParticlesIsoCone,fullName("PtSumChargedParticlesIsoCone"));
  	e.put(fTPtSumPhotonsIsoCone,fullName("PtSumPhotonsIsoCone"));
	e.put(fTDecayModeFinding,fullName("DecayModeFinding"));
  	e.put(fTVLooseIso,fullName("VLooseIso"));
  	e.put(fTLooseIso,fullName("LooseIso"));
  	e.put(fTTightIso,fullName("TightIso"));
  	e.put(fTMediumIso,fullName("MediumIso"));
/*   	e.put(fTVLooseChargedIso,fullName("VLooseChargedIso")); */
/*   	e.put(fTLooseChargedIso,fullName("LooseChargedIso")); */
/*   	e.put(fTTightChargedIso,fullName("TightChargedIso")); */
/*   	e.put(fTMediumChargedIso,fullName("MediumChargedIso")); */
/*   	e.put(fTVLooseIsoDBSumPtCorr,fullName("VLooseIsoDBSumPtCorr")); */
/*   	e.put(fTLooseIsoDBSumPtCorr,fullName("LooseIsoDBSumPtCorr")); */
/*   	e.put(fTTightIsoDBSumPtCorr,fullName("TightIsoDBSumPtCorr")); */
/*   	e.put(fTMediumIsoDBSumPtCorr,fullName("MediumIsoDBSumPtCorr")); */
/*   	e.put(fTVLooseCombinedIsoDBSumPtCorr,fullName("VLooseCombinedIsoDBSumPtCorr")); */
/*   	e.put(fTLooseCombinedIsoDBSumPtCorr,fullName("LooseCombinedIsoDBSumPtCorr")); */
/*   	e.put(fTTightCombinedIsoDBSumPtCorr,fullName("TightCombinedIsoDBSumPtCorr")); */
/*   	e.put(fTMediumCombinedIsoDBSumPtCorr,fullName("MediumCombinedIsoDBSumPtCorr")); */
  	e.put(fTLooseElectronRejection,  fullName("LooseElectronRejection"));
  	e.put(fTTightElectronRejection,  fullName("TightElectronRejection"));
  	e.put(fTMediumElectronRejection, fullName("MediumElectronRejection"));
  	e.put(fTLooseMuonRejection,      fullName("LooseMuonRejection"));
  	e.put(fTTightMuonRejection,      fullName("TightMuonRejection"));
    }else if(fType == El){
  	e.put(fTID95,fullName("ID95"));
  	e.put(fTID90,fullName("ID90"));
  	e.put(fTID85,fullName("ID85"));
  	e.put(fTID80,fullName("ID80"));
    }else if(fType == Mu){
  	e.put(fTPtErr,     fullName("PtErr"));
	e.put(fTMuNMatches,fullName("NMatches"));
    }

}

//______________________________________________________________________________
template <class LeptonType>
void LeptonFillerPat<LeptonType>::resetProducts(void) {

    fTMaxLepExc.reset(new int(0));
    fTNObjsTot .reset(new int(0));
    fTNObjs    .reset(new int(0));

    // Reset all arrays
    fTPx.reset(new std::vector<float>);
    fTPy.reset(new std::vector<float>);
    fTPz.reset(new std::vector<float>);
    fTPt.reset(new std::vector<float>);
    fTEta.reset(new std::vector<float>);
    fTPhi.reset(new std::vector<float>);
    fTE.reset(new std::vector<float>);
    fTEt.reset(new std::vector<float>);
    fTCharge.reset(new std::vector<int>);

    fTParticleIso.reset(new std::vector<float>);
    fTChargedHadronIso.reset(new std::vector<float>);
    fTNeutralHadronIso.reset(new std::vector<float>);
    fTPhotonIso.reset(new std::vector<float>);

    if(fType == Tau){
  	fTDecayMode.reset(new std::vector<int>);
	fTVz.reset(new std::vector<float>); 
	fTEmFraction.reset(new std::vector<float>); 
	fTJetPt.reset(new std::vector<float>);
	fTJetEta.reset(new std::vector<float>);
	fTJetPhi.reset(new std::vector<float>);
	fTJetMass.reset(new std::vector<float>);
	fTLeadingTkPt.reset(new std::vector<float>);
	fTLeadingNeuPt.reset(new std::vector<float>);
	fTLeadingTkHcalenergy.reset(new std::vector<float>);
	fTLeadingTkEcalenergy.reset(new std::vector<float>);
	fTNumChargedHadronsSignalCone.reset(new std::vector<int>);
	fTNumNeutralHadronsSignalCone.reset(new std::vector<int>);
	fTNumPhotonsSignalCone.reset(new std::vector<int>);
	fTNumParticlesSignalCone.reset(new std::vector<int>);
	fTNumChargedHadronsIsoCone.reset(new std::vector<int>);
	fTNumNeutralHadronsIsoCone.reset(new std::vector<int>);
	fTNumPhotonsIsolationCone.reset(new std::vector<int>);
	fTNumParticlesIsolationCone.reset(new std::vector<int>);
	fTPtSumChargedParticlesIsoCone.reset(new std::vector<float>);
	fTPtSumPhotonsIsoCone.reset(new std::vector<float>);
	fTDecayModeFinding.reset(new std::vector<float>);
	fTVLooseIso.reset(new std::vector<float>);
	fTLooseIso.reset(new std::vector<float>);
	fTTightIso.reset(new std::vector<float>);
	fTMediumIso.reset(new std::vector<float>);
/* 	fTVLooseChargedIso.reset(new std::vector<float>); */
/* 	fTLooseChargedIso.reset(new std::vector<float>); */
/* 	fTTightChargedIso.reset(new std::vector<float>); */
/* 	fTMediumChargedIso.reset(new std::vector<float>); */
/* 	fTVLooseIsoDBSumPtCorr.reset(new std::vector<float>); */
/* 	fTLooseIsoDBSumPtCorr.reset(new std::vector<float>); */
/* 	fTTightIsoDBSumPtCorr.reset(new std::vector<float>); */
/* 	fTMediumIsoDBSumPtCorr.reset(new std::vector<float>); */
/* 	fTVLooseCombinedIsoDBSumPtCorr.reset(new std::vector<float>); */
/* 	fTLooseCombinedIsoDBSumPtCorr.reset(new std::vector<float>); */
/* 	fTTightCombinedIsoDBSumPtCorr.reset(new std::vector<float>); */
/* 	fTMediumCombinedIsoDBSumPtCorr.reset(new std::vector<float>); */
	fTLooseElectronRejection.reset(new std::vector<float>);
	fTTightElectronRejection.reset(new std::vector<float>);
	fTMediumElectronRejection.reset(new std::vector<float>);
	fTLooseMuonRejection.reset(new std::vector<float>);
	fTTightMuonRejection.reset(new std::vector<float>);
    }else if(fType == El){
  	fTID95.reset(new std::vector<int>);
  	fTID90.reset(new std::vector<int>);
  	fTID85.reset(new std::vector<int>);
  	fTID80.reset(new std::vector<int>);
    }else if(fType == Mu){
  	fTMuNMatches.reset(new std::vector<int>);
	fTPtErr.reset(new std::vector<float>);
    }
  

}

//________________________________________________________________________________________
template <>
void LeptonFillerPat<pat::Tau>::getSpecific(const pat::Tau& lepton){
    // speficic for PFTaus
    fTDecayMode  ->push_back( lepton.decayMode() );	
    fTVz         ->push_back( lepton.vz() );  
    fTEmFraction ->push_back( lepton.emFraction() ); 
    fTJetPt      ->push_back( lepton.pfJetRef().get()->pt() );
    fTJetEta     ->push_back( lepton.pfJetRef().get()->eta() );
    fTJetPhi     ->push_back( lepton.pfJetRef().get()->phi() );
    fTJetMass    ->push_back( lepton.pfJetRef().get()->mass() );

    fTLeadingTkPt        ->push_back( (lepton.leadPFChargedHadrCand())->pt() );
    fTLeadingNeuPt       ->push_back( (lepton.leadPFNeutralCand().isNonnull() ? lepton.leadPFNeutralCand()->pt() : 0.) );
    fTLeadingTkHcalenergy->push_back( lepton.leadPFChargedHadrCand()->hcalEnergy() );
    fTLeadingTkEcalenergy->push_back( lepton.leadPFChargedHadrCand()->ecalEnergy() );

    fTNumChargedHadronsSignalCone->push_back( lepton.signalPFChargedHadrCands().size() );
    fTNumNeutralHadronsSignalCone->push_back( lepton.signalPFNeutrHadrCands().size() );
    fTNumPhotonsSignalCone->push_back( lepton.signalPFGammaCands().size() );
    fTNumParticlesSignalCone->push_back( lepton.signalPFCands().size() );

    fTNumChargedHadronsIsoCone->push_back( lepton.isolationPFChargedHadrCands().size() );
    fTNumNeutralHadronsIsoCone->push_back( lepton.isolationPFNeutrHadrCands().size() );
    fTNumPhotonsIsolationCone->push_back( lepton.isolationPFGammaCands().size() );
    fTNumParticlesIsolationCone->push_back( lepton.isolationPFCands().size() );
    fTPtSumChargedParticlesIsoCone->push_back( lepton.isolationPFChargedHadrCandsPtSum() );
    fTPtSumPhotonsIsoCone->push_back( lepton.isolationPFGammaCandsEtSum() );

    fTDecayModeFinding->push_back( lepton.tauID("decayModeFinding") );
    fTVLooseIso       ->push_back( lepton.tauID("byVLooseIsolation") );
    fTLooseIso        ->push_back( lepton.tauID("byLooseIsolation") );
    fTTightIso        ->push_back( lepton.tauID("byTightIsolation") );
    fTMediumIso       ->push_back( lepton.tauID("byMediumIsolation") );
/* 	fTVLooseChargedIso->push_back( lepton.tauID("byVLooseChargedIsolation") ); */
/* 	fTLooseChargedIso->push_back( lepton.tauID("byLooseChargedIsolation") ); */
/* 	fTTightChargedIso->push_back( lepton.tauID("byTightChargedIsolation") ); */
/* 	fTMediumChargedIso->push_back( lepton.tauID("byMediumChargedIsolation") ); */
/*         fTVLooseIsoDBSumPtCorr->push_back( lepton.tauID("byVLooseIsolationDeltaBetaCorr") ); */
/* 	fTLooseIsoDBSumPtCorr->push_back( lepton.tauID("byLooseIsolationDeltaBetaCorr") ); */
/* 	fTTightIsoDBSumPtCorr->push_back( lepton.tauID("byTightIsolationDeltaBetaCorr") ); */
/* 	fTMediumIsoDBSumPtCorr->push_back( lepton.tauID("byMediumIsolationDeltaBetaCorr") ); */
/* 	fTVLooseCombinedIsoDBSumPtCorr->push_back( lepton.tauID("byVLooseCombinedIsolationDeltaBetaCorr") ); */
/* 	fTLooseCombinedIsoDBSumPtCorr->push_back( lepton.tauID("byLooseCombinedIsolationDeltaBetaCorr") ); */
/* 	fTTightCombinedIsoDBSumPtCorr->push_back( lepton.tauID("byTightCombinedIsolationDeltaBetaCorr") ); */
/* 	fTMediumCombinedIsoDBSumPtCorr->push_back( lepton.tauID("byMediumCombinedIsolationDeltaBetaCorr") ); */

    fTLooseElectronRejection ->push_back( lepton.tauID("againstElectronLoose") );
    fTTightElectronRejection ->push_back( lepton.tauID("againstElectronTight") );
    fTMediumElectronRejection->push_back( lepton.tauID("againstElectronMedium") );
    
    fTLooseMuonRejection->push_back( lepton.tauID("againstMuonLoose") );
    fTTightMuonRejection->push_back( lepton.tauID("againstMuonTight") );
    return;
}

template <>
void LeptonFillerPat<pat::Electron>::getSpecific(const pat::Electron& lepton){
    // speficic for PFElectrons
    fTID95->push_back( lepton.electronID("simpleEleId95cIso") );
    fTID90->push_back( lepton.electronID("simpleEleId90cIso") );
    fTID85->push_back( lepton.electronID("simpleEleId85cIso") );
    fTID80->push_back( lepton.electronID("simpleEleId80cIso") );
    return;
}

template <>
void LeptonFillerPat<pat::Muon>::getSpecific(const pat::Muon& lepton){
    // speficic for PFMuon
    fTMuNMatches->push_back( lepton.numberOfMatches() );
    fTPtErr     ->push_back( lepton.globalTrack()->ptError() );
    return;
}

#endif
