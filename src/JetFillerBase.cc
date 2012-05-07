#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/JPTJetCollection.h"
#include "DataFormats/JetReco/interface/JetID.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "DataFormats/JetReco/interface/JetTracksAssociation.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "DiLeptonAnalysis/NTupleProducer/interface/JetFillerBase.h"


//________________________________________________________________________________________
JetFillerBase::JetFillerBase( const edm::ParameterSet& cfg, const bool& isRealData )
    : FillerBase(cfg,isRealData)
{
	
    // Retrieve configuration parameters
    std::string jettype = cfg.getParameter<edm::InputTag>("tag").label();

    // parse fTag label for jet-type
    if      (std::string::npos != jettype.find("PF"))   setJetType(PF);
    else if (std::string::npos != jettype.find("Calo")) setJetType(CALO);
    else if (std::string::npos != jettype.find("JPT"))  setJetType(JPT);  
    else {
        setJetType(unknown);
        edm::LogWarning("NTP") << "!! Don't know JetType !!" << jettype;
    }

    //FIXME: could be set from configuration file...
    gMaxnobjs = 100;

}

//________________________________________________________________________________________
const std::vector<filler::PPair> JetFillerBase::declareProducts(void) {

    addProduct("NJets",   typeid(*fTNObjs));
    addProduct("JPx",     typeid(*fTPx));
    addProduct("JPy",     typeid(*fTPy));
    addProduct("JPz",     typeid(*fTPz));
    addProduct("JPt",     typeid(*fTPt));
    addProduct("JE",      typeid(*fTE));
    addProduct("JEt",     typeid(*fTEt));
    addProduct("JEta",    typeid(*fTEta));
    addProduct("JPhi",    typeid(*fTPhi));
    addProduct("JFlavour",typeid(*fTFlavour));

    addProduct("JScale",  typeid(*fTScale));
    addProduct("JL1FastJetScale", typeid(*fTL1FastJetScale));
    addProduct("JArea",   typeid(*fTArea));
    addProduct("JIDLoose",typeid(*fTIDLoose));
    for ( std::vector<std::string>::const_iterator it = fBtagNames.begin();
        it != fBtagNames.end(); ++it ) {
        addProduct(("J"+(*it)).c_str(), typeid((*fTJbTagProb[0])));

    }
    if(fJetType==CALO) {
        addProduct("JNConstituents", typeid(*fTNConstituents));
        addProduct("JNAssoTracks",   typeid(*fTNAssoTracks));
        addProduct("JChfrac",        typeid(*fTChfrac));
        addProduct("JEMfrac",        typeid(*fTEMfrac));
        addProduct("JIDHPD",         typeid(*fTID_HPD));
        addProduct("JIDRBX",         typeid(*fTID_RBX));
        addProduct("JIDn90Hits",     typeid(*fTID_n90Hits));
        addProduct("Jn90",           typeid(*fTn90));
        addProduct("JIDresEMF",      typeid(*fTID_resEMF));
    } else if (fJetType==JPT) {
        addProduct("JChMult",    typeid(*fTChMult));
        addProduct("JIDHPD",    typeid(*fTID_HPD));
        addProduct("JIDRBX",    typeid(*fTID_RBX));
        addProduct("JIDn90Hits",typeid(*fTID_n90Hits));
        addProduct("JIDresEMF", typeid(*fTID_resEMF));
    } else if (fJetType==PF) {
        addProduct("JNConstituents", typeid(*fTNConstituents));
        addProduct("JChMult",    typeid(*fTChMult));
        addProduct("JNeuMult",   typeid(*fTNeuMult));
        addProduct("JChHadfrac", typeid(*fTChHadfrac));
        addProduct("JNeuHadfrac",typeid(*fTNeuHadfrac));
        addProduct("JChEmfrac",  typeid(*fTChEmfrac));
        addProduct("JNeuEmfrac", typeid(*fTNeuEmfrac));
        addProduct("JChMufrac",  typeid(*fTChMufrac));
        addProduct("JPhofrac",   typeid(*fTPhofrac));
        addProduct("JHFHadfrac", typeid(*fTHFHadfrac));
        addProduct("JHFEMfrac",  typeid(*fTHFEMfrac));
    }

    return typeList;
}

//________________________________________________________________________________________
void JetFillerBase::resetProducts(void) {

    fTNObjs.reset(new int(0));
	
    fTPx     .reset(new std::vector<float>);
    fTPy     .reset(new std::vector<float>);
    fTPz     .reset(new std::vector<float>);
    fTPt     .reset(new std::vector<float>);
    fTE      .reset(new std::vector<float>);
    fTEt     .reset(new std::vector<float>);
    fTEta    .reset(new std::vector<float>);
    fTPhi    .reset(new std::vector<float>);
    fTArea   .reset(new std::vector<float>);
    fTScale  .reset(new std::vector<float>);
    fTL1FastJetScale.reset(new std::vector<float>);
    fTFlavour.reset(new std::vector<int>);
    size_t ibtag = 0;
    for ( std::vector<std::string>::const_iterator it = fBtagNames.begin();
        it != fBtagNames.end(); ++it ) {
        fTJbTagProb[ibtag++].reset(new std::vector<float> );
    }
    fTIDLoose.reset(new std::vector<int>);
    if (fJetType==CALO) {	
        fTNConstituents.reset(new std::vector<int>);		
        fTNAssoTracks .reset(new std::vector<int>);
        fTChfrac      .reset(new std::vector<float>);
        fTEMfrac       .reset(new std::vector<float>);
        fTID_RBX       .reset(new std::vector<float>);
        fTID_HPD       .reset(new std::vector<float>);
        fTID_n90Hits   .reset(new std::vector<float>);
        fTn90          .reset(new std::vector<int>);
        fTID_resEMF    .reset(new std::vector<float>);
    } else if (fJetType==JPT) {	
        fTID_RBX       .reset(new std::vector<float>);
        fTID_HPD       .reset(new std::vector<float>);
        fTID_n90Hits   .reset(new std::vector<float>);
        fTID_resEMF    .reset(new std::vector<float>);
        fTChMult       .reset(new std::vector<int>);
    } else if (fJetType==PF) {
        fTNConstituents.reset(new std::vector<int>);
        fTChMult       .reset(new std::vector<int>);
        fTNeuMult      .reset(new std::vector<int>);
        fTChHadfrac    .reset(new std::vector<float>);
        fTNeuHadfrac   .reset(new std::vector<float>);
        fTChEmfrac     .reset(new std::vector<float>);
        fTNeuEmfrac    .reset(new std::vector<float>);
        fTChMufrac     .reset(new std::vector<float>);
        fTPhofrac      .reset(new std::vector<float>);
        fTHFHadfrac    .reset(new std::vector<float>);
        fTHFEMfrac     .reset(new std::vector<float>);
    }
}


//______________________________________________________________________________
void JetFillerBase::putProducts( edm::Event& e ) { 

    e.put(fTNObjs,fullName("NJets"));
	
    e.put(fTPx, fullName("JPx"));
    e.put(fTPy, fullName("JPy"));
    e.put(fTPz, fullName("JPz"));
    e.put(fTPt, fullName("JPt"));
    e.put(fTE,  fullName("JE"));
    e.put(fTEt, fullName("JEt"));
    e.put(fTEta,fullName("JEta"));
    e.put(fTPhi,fullName("JPhi"));
    e.put(fTFlavour,fullName("JFlavour"));
    e.put(fTScale,  fullName("JScale"));
    e.put(fTL1FastJetScale,fullName("JL1FastJetScale"));
    e.put(fTArea,   fullName("JArea"));
    e.put(fTIDLoose,fullName("JIDLoose"));
    size_t ibtag = 0;
    for ( std::vector<std::string>::const_iterator it = fBtagNames.begin();
        it != fBtagNames.end(); ++it ) {
        e.put(fTJbTagProb[ibtag++], fullName(("J"+(*it)).c_str()));
    }
    if (fJetType==CALO) {	
        e.put(fTNConstituents,fullName("JNConstituents")); 
        e.put(fTNAssoTracks,  fullName("JNAssoTracks"));
        e.put(fTChfrac,       fullName("JChfrac"));
        e.put(fTEMfrac,       fullName("JEMfrac"));
        e.put(fTID_HPD,       fullName("JIDHPD"));
        e.put(fTID_RBX,       fullName("JIDRBX"));
        e.put(fTID_n90Hits,   fullName("JIDn90Hits"));
        e.put(fTn90,          fullName("Jn90"));
        e.put(fTID_resEMF,    fullName("JIDresEMF"));
    } else if (fJetType==JPT) {	
        e.put(fTChMult,       fullName("JChMult"));
        e.put(fTID_HPD,       fullName("JIDHPD"));
        e.put(fTID_RBX,       fullName("JIDRBX"));
        e.put(fTID_n90Hits,   fullName("JIDn90Hits"));
        e.put(fTID_resEMF,    fullName("JIDresEMF"));
    } else if (fJetType==PF) {
        e.put(fTNConstituents,fullName("JNConstituents"));
        e.put(fTChMult,       fullName("JChMult"));
        e.put(fTNeuMult,      fullName("JNeuMult"));
        e.put(fTChHadfrac,    fullName("JChHadfrac"));
        e.put(fTNeuHadfrac,   fullName("JNeuHadfrac"));
        e.put(fTChEmfrac,     fullName("JChEmfrac"));
        e.put(fTNeuEmfrac,    fullName("JNeuEmfrac"));
        e.put(fTChMufrac,     fullName("JChMufrac"));
        e.put(fTPhofrac,      fullName("JPhofrac"));
        e.put(fTHFHadfrac,    fullName("JHFHadfrac"));
        e.put(fTHFEMfrac,     fullName("JHFEMfrac"));
    }
}
