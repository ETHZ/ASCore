#include "helper/Utilities.hh"
#include "MassAnalysis.hh"
#include "TLorentzVector.h"
#include <sstream>

using namespace std;

MassAnalysis::MassAnalysis(TreeReader *tr) : MultiplicityAnalysisBase(tr){
	Util::SetStyle();	
}

MassAnalysis::~MassAnalysis(){
}

void MassAnalysis::Begin(const char* filename){
	// Define the output file of histograms
	fHistFile = new TFile(fOutputDir + TString(filename), "RECREATE");

	// book tree
	fMT2tree = new MT2tree();
	BookTree();
}

void MassAnalysis::Analyze(){	

	// ---------------------------------------------------
	// Initialize fElecs, fJetsLoose, fBJets, fMuons, fLeptConfig 
	InitializeEvent();
	// ----------------------------------------------------

	// jet dataset
	double run[7]   = {147216, 147217, 147390, 147929, 148031, 148032, 147048};
	double lumi[7]  = {47, 179, 279, 378, 771, 144, 112};
	double event[7] = {35190778, 131464723, 242410506, 352665576, 599244226, 106359451,  72316051};

	// multijet
	//double run[21]   = {148822, 148829, 148862, 148862, 148864, 148864, 148953, 149011, 149011, 149011, 149181, 
	//                    149181, 149181, 149181, 149291, 149291, 148864, 148864, 148864, 149058, 149181};
	//double lumi[21]  = {440, 244, 289, 497, 460, 514,  49,  67, 424, 461, 257, 305, 707, 941, 436, 539,  10,  47, 275,  64, 985};
	//double event[21] = {454994856, 239439139, 434297633, 729330524, 529426262, 586943934,  70626194, 105062641, 
	//                    619833331, 668579778,  87276054, 157837774, 670139805, 925116042, 466605433, 569274059,   9986520,  56130889, 325562876,  76943429, 970118251};
	
	// LM1
	//double event[50] = {84407, 84410, 84411, 84413, 84415, 84429, 84431, 84434, 84437, 84442, 84443, 84444, 84445, 
	//                    84447, 84452, 84455, 84458, 84460, 84461, 84464, 84465, 84470, 84472, 84484, 84497, 84500, 
	//                    84504, 84523, 84524, 84527, 84531, 84534, 84547, 84553, 84562, 84568, 84576, 84577, 84579, 
	//                    84591, 84594, 84595, 84604, 84605, 84614, 84623, 84625, 84629, 84633, 84635};
	

//	for(int i=0; i<7; ++i){
//		if(fTR->Run == run[i] && fTR->LumiSection==lumi[i] && fTR->Event==event[i]){
//			UserAnalysisBase::EventPrint();
//		}
//	}

	// --------------------------------------------------------------------
	// check if event passes selection cuts and trigger requirements
	if(! IsSelectedEvent()){return;}
	// --------------------------------------------------------------------
		
	// reset tree variables
	ResetTree();
	
	// -----------------------------------------------------------
	// fill fHemiObjects1,2,3
	for(int i=0; i<gNHemispheres; ++i){ fHemiObjects[i].Reset();}
	   // ************************************************* //
	   // changes here affect misc.MT2 type variables !!!!  //
	   // ************************************************  //
	   // seed: max inv mass, assoc: minimal lund distance, no dR cut (i.e. third parameter is 0)
	GetMT2Variables(2, 3, 0.0, 20, 2.4, fHemiObjects[0]); 
	
     	// MT2 for two leading jets, maxdR=0.1, i.e. smaller than jet cone size
	GetMT2Variables(4, 3, 0.1, 20, 2.4, fHemiObjects[1]); 

	// seed: leading jets, assoc: minimal lund distance, maxdR=1
	GetMT2Variables(4, 3, 1.0, 20, 2.4, fHemiObjects[2]); 

	// minimizing deltaHT
	bool minimizeDHT(true);
	GetMT2Variables(minimizeDHT, 20, 2.4, fHemiObjects[3]); 
	// -----------------------------------------------------------


	// Fill Tree !! has to be called at the very end !!
	MassAnalysis::FillTree();
}

// ***********************************************************************************************
void MassAnalysis::BookTree(){
	           
	fATree = new TTree("MassTree", "MassTree");
	fATree->Branch("MT2tree" , "MT2tree" , &fMT2tree);

}

void MassAnalysis::ResetTree(){
        fMT2tree->Reset();
}

void MassAnalysis::FillTree(){
	// ------------------------------------------------------------------
	// fill misc 
	fMT2tree->misc.Run                 = fTR->Run;
	fMT2tree->misc.Event		   = fTR->Event;
	fMT2tree->misc.LumiSection	   = fTR->LumiSection;
	fMT2tree->misc.LeptConfig          = (int) fLeptConfig;
	fMT2tree->misc.HBHENoiseFlag	   = fTR->HBHENoiseFlag;
	fMT2tree->misc.Vectorsumpt	   = fVectorSumPt;
	fMT2tree->misc.HT                  = fHT;
	fMT2tree->misc.PFMETsign	   = (fTR->PFMET)/sqrt(fTR->SumEt);
	fMT2tree->misc.EcalDeadCellBEFlag  = fTR->EcalDeadCellBEFlag;
	fMT2tree->misc.NECALGapClusters    = fTR->NECALGapClusters;
	for(int i=0; i<fMT2tree->misc.NECALGapClusters; ++i){
		fMT2tree->misc.EcalGapClusterSize[i] = fTR->EcalGapClusterSize[i];
		fMT2tree->misc.EcalGapBE[i]          = fTR->EcalGapBE[i];
	}
	fMT2tree->misc.MT2                 = fHemiObjects[0].MT2;    // note: this is a bit dangerous, 
	fMT2tree->misc.MCT                 = fHemiObjects[0].MCT;
	fMT2tree->misc.MT2leading          = fHemiObjects[1].MT2;    //       as the definition of fHemiObjects1,2,3,4
	fMT2tree->misc.MT2noISR            = fHemiObjects[2].MT2;    //       is interchangable
	fMT2tree->misc.AlphaT              = fHemiObjects[3].alphaT; 

	fMT2tree->misc.MET                 = fTR->PFMET;
	fMT2tree->misc.METPhi              = fTR->PFMETphi;

	// NVertices
	int nvertex=0;
	for(int i=0; i<fTR->NVrtx; ++i){
		if(fabs(fTR->VrtxZ[i]) > 24) continue;
		if(sqrt( (fTR->VrtxX[i])*(fTR->VrtxX[i]) + (fTR->VrtxY[i])*(fTR->VrtxY[i])) > 2) continue;
		if(fTR->VrtxNdof[i]<=4) continue;
		nvertex++;
	}
	fMT2tree->misc.NVertices=nvertex;

	// get number of jets with pt > 20 and |eta| < 5
	int npfjets=0;
	for(int i=0; i<fTR->PFNJets; ++i){
		if(fTR->PFJPt[i]>20 && fabs(fTR->PFJEta[i])<5.0){
			npfjets++;	
		}
	}	
	fMT2tree->misc.NJetsEta5Pt20 = npfjets;

	// ----------------------------------------------------------------
	// fill MT2Hemi
	for(int h=0; h<gNHemispheres; ++h){
		fMT2tree->hemi[h].assoc_method  = fHemiObjects[h].assoc;
		fMT2tree->hemi[h].seed_method   = fHemiObjects[h].seed;
		fMT2tree->hemi[h].MT2           = fHemiObjects[h].MT2;
		fMT2tree->hemi[h].MCT           = fHemiObjects[h].MCT;
		fMT2tree->hemi[h].AlphaT        = fHemiObjects[h].alphaT;
		fMT2tree->hemi[h].minDHT        = fHemiObjects[h].minDHT;
		fMT2tree->hemi[h].lv1           = fHemiObjects[h].pjet1;
		fMT2tree->hemi[h].lv2           = fHemiObjects[h].pjet2;
		fMT2tree->hemi[h].UTM           = fHemiObjects[h].UTM;

		int jcount1=0, jcount2=0, elecount1=0, elecount2=0, muocount1=0, muocount2=0; 
		for(int i=0; i<fHemiObjects[h].objects.size(); ++i){
			HemiObject obj=fHemiObjects[h].objects[i];
			if(obj.type=="jet" && obj.hemi==1){
				fMT2tree->hemi[h].jindices1[jcount1]=obj.index;
				jcount1++;
			} else if(obj.type=="jet" && obj.hemi==2){
				fMT2tree->hemi[h].jindices2[jcount2]=obj.index;
				jcount2++;
			} else if(obj.type=="ele" && obj.hemi==1){
				fMT2tree->hemi[h].eleindices1[elecount1]=obj.index;
				elecount1++;
			} else if(obj.type=="ele" && obj.hemi==2){
				fMT2tree->hemi[h].eleindices2[elecount2]=obj.index;
				elecount2++;
			} else if(obj.type=="muo" && obj.hemi==1){
				fMT2tree->hemi[h].muoindices1[muocount1]=obj.index;
				muocount1++;
			} else if(obj.type=="muo" && obj.hemi==2){
				fMT2tree->hemi[h].muoindices2[muocount2]=obj.index;
				muocount2++;
			}
		}
	}
	
	// ---------------------------------------------------------------
	// Set NJets, NElecs, NMuons
	fMT2tree->SetNJets         ((Int_t)fJets.size());
	fMT2tree->SetNJetsIDLoose  ((Int_t)fJetsLoose.size());
	fMT2tree->SetNJetsIDMedium ((Int_t)fJetsMedium.size());
	fMT2tree->SetNJetsIDTight  ((Int_t)fJetsTight.size());
	fMT2tree->SetNEles         ((Int_t)fElecs.size());
	fMT2tree->SetNElesLoose    ((Int_t)fElecsLoose.size());
	fMT2tree->SetNMuons        ((Int_t)fMuons.size());
	fMT2tree->SetNMuonsLoose   ((Int_t)fMuonsLoose.size());
	
	// ---------------------------------------------------------------
	// Fill jets 4-momenta & ID's
	for(int i=0; i<fJets.size(); ++i) {
	  	fMT2tree->jet[i].lv.SetPxPyPzE( fTR->PFJPx[fJets[i]], fTR->PFJPy[fJets[i]], 
					fTR->PFJPz[fJets[i]], fTR->PFJE [fJets[i]]); //SetLV(GetJet4Momenta(fJets[i]));
		// b-tag info now should be available
	 	fMT2tree->jet[i].bTagProbTCHE  =  fTR->PFJbTagProbTkCntHighEff [fJets[i]];
		fMT2tree->jet[i].bTagProbTCHP  =  fTR->PFJbTagProbTkCntHighPur [fJets[i]];
		fMT2tree->jet[i].bTagProbSSVHE =  fTR->PFJbTagProbSimpSVHighEff[fJets[i]];
		fMT2tree->jet[i].bTagProbSSVHP =  fTR->PFJbTagProbSimpSVHighPur[fJets[i]];
	  
		// Jet id variables
		if(IsGoodBasicPFJet (fJets[i], 20., 2.4))  fMT2tree->jet[i].isPFIDLoose =true;
		if(IsGoodPFJetMedium(fJets[i], 20., 2.4))  fMT2tree->jet[i].isPFIDMedium=true;
		if(IsGoodPFJetTight (fJets[i], 20., 2.4))  fMT2tree->jet[i].isPFIDTight =true;
		fMT2tree->jet[i].ChHadFrac      = fTR->PFJChHadfrac     [fJets[i]];	
		fMT2tree->jet[i].NeuHadFrac     = fTR->PFJNeuHadfrac    [fJets[i]];
		fMT2tree->jet[i].ChEmFrac       = fTR->PFJChEmfrac      [fJets[i]];
		fMT2tree->jet[i].NeuEmFrac      = fTR->PFJNeuEmfrac     [fJets[i]];
		fMT2tree->jet[i].ChMult         = fTR->PFJChMult        [fJets[i]];
		fMT2tree->jet[i].NeuMult        = fTR->PFJNeuMult       [fJets[i]];
		fMT2tree->jet[i].NConstituents  = fTR->PFJNConstituents [fJets[i]];
	}
	
	// -----------------------------------------------------------------
	// Fill leptons 4-momenta & tight_flag
	for(int i=0; i<fElecsLoose.size(); ++i) {
	  	fMT2tree->ele[i].lv.SetPtEtaPhiE(fTR->ElPt [fElecsLoose[i]], fTR->ElEta[fElecsLoose[i]], 
				      fTR->ElPhi[fElecsLoose[i]], fTR->ElE  [fElecsLoose[i]]); // = GetEle4Momenta(fElecs[i]);
	  	for(int l=0; l<fElecs.size(); ++l) {
	  		if(fElecs[l]==fElecsLoose[i]) fMT2tree->ele[i].isTight=true;
	  	}
	}
	for(int i=0; i<fMuonsLoose.size(); ++i) {
	  	fMT2tree->muo[i].lv.SetPtEtaPhiM(fTR->MuPt [fMuonsLoose[i]], fTR->MuEta[fMuonsLoose[i]], 
				      fTR->MuPhi[fMuonsLoose[i]], 0.106);                     // = GetMuo4Momenta(fMuons[i]);
	  	for(int l=0; l<fMuons.size(); ++l) {
	  		if(fMuons[l]==fMuonsLoose[i]) fMT2tree->muo[i].isTight=true;
	  	}
	}
	
	// -------------------------------------------------------------------
	// Genleptons
	for(int i=0; i<fTR->NGenLeptons; ++i){
		double mass=0;
		if     (abs(fTR->GenLeptonID[i]) == 15) mass=1.776;
		else if(abs(fTR->GenLeptonID[i]) == 13) mass=0.106;

		fMT2tree->genlept[i].lv.SetPtEtaPhiM(fTR->GenLeptonPt[i], fTR->GenLeptonEta[i], fTR->GenLeptonPhi[i], mass);
		fMT2tree->genlept[i].ID       = fTR->GenLeptonID[i];
		fMT2tree->genlept[i].MID      = fTR->GenLeptonMID[i];
		fMT2tree->genlept[i].MStatus  = fTR->GenLeptonMStatus[i];
		fMT2tree->genlept[i].GMID     = fTR->GenLeptonGMID[i];
		fMT2tree->genlept[i].GMStatus = fTR->GenLeptonGMStatus[i];
	}


	// --------------------------------------------------------------------
	// MET, MHT and MPT
	fMT2tree->pfmet[0].SetPtEtaPhiM(fTR->PFMET,0,fTR->PFMETphi,0);

	// Fill vector sum of tracks
	TVector3 tracks(0.,0.,0.);
	for(int i=0; i< fTR->NTracks; ++i){
		TVector3 track;
		track.SetPtEtaPhi(fabs(fTR->TrkPt[i]), fTR->TrkEta[i], fTR->TrkPhi[i]);
		tracks += track;
	}	
	fMT2tree->MPT[0].SetXYZM(-tracks.Px(), -tracks.Py(), 0, 0);

	// DPhiMhtMpt: should be replaced by method to calculate it on the fly. 
	fMT2tree->misc.DPhiMhtMpt=Util::DeltaPhi(fMHTphi, fMT2tree->MPT[0].Phi());	
	
	// fill MHT
	TVector3 MHTidloose(0., 0., 0.);
	TVector3 MHTall(0., 0., 0.);
        for(int i=0; i<fTR->PFNJets; ++i) {
		TVector3 jet;
		jet.SetPtEtaPhi(fTR->PFJPt[i], fTR->PFJEta[i], fTR->PFJPhi[i]);
		MHTall += jet;
		if(! IsGoodBasicPFJet (i, 20., 2.4)) continue;
		MHTidloose += jet;
	}	
	fMT2tree->MHTloose[0].SetXYZM(-MHTidloose.Px(), -MHTidloose.Py(), 0, 0);
	fMT2tree->MHT     [0].SetXYZM(-MHTall.Px()    , -MHTall.Py()    , 0, 0);


	// -----------------------------------------------------------------------
	// fill tree
	fATree            ->Fill();
}

// ************************************************************************************************************************************* 
void MassAnalysis::GetMT2Variables(int hemi_seed, int hemi_assoc, double maxDR, double minJPt, double maxJEta, HemiObjects& hemiobject){

	// clear hemiobject
	hemiobject.objects .clear();
	hemiobject.pjet1   .Clear();
	hemiobject.pjet2   .Clear();
	hemiobject.UTM     .Clear();
	hemiobject.MT2     =-999.99;
	hemiobject.MCT     =-999.99;

	// fill stuff
	hemiobject.assoc   = hemi_assoc;
	hemiobject.seed    = hemi_seed;
	hemiobject.minDHT  =-999.99;
	hemiobject.alphaT  =-999.99;


	// make pseudojets with hemispheres
	vector<float> px, py, pz, E;
	for(int i=0; i<fJets.size(); ++i){
		if(IsGoodBasicPFJet(i, minJPt, maxJEta)==false) continue;
		px.push_back(fTR->PFJPx[fJets[i]]);
		py.push_back(fTR->PFJPy[fJets[i]]);
		pz.push_back(fTR->PFJPz[fJets[i]]);
		 E.push_back(fTR->PFJE[fJets[i]]);
		 fHemiObject.type="jet"; fHemiObject.index=fJets[i]; fHemiObject.hemi=0;
		 hemiobject.objects.push_back(fHemiObject);
	}
	for(int i=0; i< fElecs.size(); ++i){
		px.push_back(fTR->ElPx[fElecs[i]]);
		py.push_back(fTR->ElPy[fElecs[i]]);
		pz.push_back(fTR->ElPz[fElecs[i]]);
		 E.push_back(fTR->ElE [fElecs[i]]);
		 fHemiObject.type="ele"; fHemiObject.index=fElecs[i]; fHemiObject.hemi=0;
		 hemiobject.objects.push_back(fHemiObject);
	}
	for(int i=0; i< fMuons.size(); ++i){
		px.push_back(fTR->MuPx[fMuons[i]]);
		py.push_back(fTR->MuPy[fMuons[i]]);
		pz.push_back(fTR->MuPz[fMuons[i]]);
		 E.push_back(fTR->MuE [fMuons[i]]);
		 fHemiObject.type="muo"; fHemiObject.index=fMuons[i]; fHemiObject.hemi=0;
		 hemiobject.objects.push_back(fHemiObject);
	}
	// get hemispheres 
	Hemisphere* hemi      = new Hemisphere(px, py, pz, E, hemi_seed, hemi_assoc);
	if(maxDR>0) hemi      -> RejectISRDRmax(maxDR);  // consider only jets withing dR < maxDR from hemi-axix for pseudojets
	hemi                  -> SetDebug(0);

	vector<int> grouping  = hemi ->getGrouping();
	
	TLorentzVector pseudojet1(0.,0.,0.,0.);
	TLorentzVector pseudojet2(0.,0.,0.,0.);
	
	for(int i=0; i<px.size(); ++i){
		if(grouping[i]==1){
			pseudojet1.SetPx(pseudojet1.Px() + px[i]);
			pseudojet1.SetPy(pseudojet1.Py() + py[i]);
			pseudojet1.SetPz(pseudojet1.Pz() + pz[i]);
			pseudojet1.SetE( pseudojet1.E()  + E[i]);
			hemiobject.objects[i].hemi=1;
		}else if(grouping[i] == 2){
			pseudojet2.SetPx(pseudojet2.Px() + px[i]);
			pseudojet2.SetPy(pseudojet2.Py() + py[i]);
			pseudojet2.SetPz(pseudojet2.Pz() + pz[i]);
			pseudojet2.SetE( pseudojet2.E()  + E[i]);
			hemiobject.objects[i].hemi=2;
		}
	}
	delete hemi;

	TLorentzVector pmiss;
	pmiss.SetPx(fTR->PFMETpx);
	pmiss.SetPy(fTR->PFMETpy);

	// fill UTM
	hemiobject.UTM          = - pmiss - pseudojet1 - pseudojet2;

	// fill MT2
	hemiobject.MT2   =GetMT2(pseudojet1, 0., pseudojet2, 0., pmiss, 0.);
	
	// fill MCT
	TVector2 pmiss_vector2;
	pmiss_vector2.Set(fTR->PFMETpx, fTR->PFMETpy);
	TLorentzVector downstream(0.,0.,0.,0.); // no particles are downstream, i.e. not selected jets are upstream. 
	hemiobject.MCT   =GetMCTcorr(pseudojet1, pseudojet2, downstream, pmiss_vector2);
	
	// fill pseudojets
	hemiobject.pjet1 =pseudojet1;
	hemiobject.pjet2 =pseudojet2;

	// --------------------------------
	// testing: comparing to ToveyMT2	
	fMCT = new TMctLib;
	double MT2Tovey  =fMCT->mt2(pseudojet1, pseudojet2, downstream, pmiss_vector2, 14000, 0);
	delete fMCT;
	if(fVerbose > 1) cout << "Event " << fTR->Event <<" MT2 " << hemiobject.MT2 
	                      << " MT2Tovey " << MT2Tovey << endl;
	// --------------------------------
}

void MassAnalysis::GetMT2Variables(bool minimizeDHT,  double minJPt, double maxJEta, HemiObjects& hemiobject){

	// FIXME: save info about which jet is in which hemisphere e.g. fHemiObject.hemi 
	if(minimizeDHT !=1){ cout << "MassAnalysis::GetMT2Variables: ERROR: got minimizeDHT!=1" << endl; return; }

	// clear hemiobject
	hemiobject.objects .clear();
	hemiobject.pjet1   .Clear();
	hemiobject.pjet2   .Clear();
	hemiobject.UTM     .Clear();
	hemiobject.MT2     =-999.99;
	hemiobject.MCT     =-999.99;

	// fill stuff
	hemiobject.assoc   = -1;
	hemiobject.seed    = -1;


	// make pseudojets with hemispheres
	vector<TLorentzVector> p4s;
	for(int i=0; i<fJets.size(); ++i){
		if(IsGoodBasicPFJet(i, minJPt, maxJEta)==false) continue;
		TLorentzVector v; 
		v.SetPxPyPzE(fTR->PFJPx[fJets[i]], fTR->PFJPy[fJets[i]], fTR->PFJPz[fJets[i]], fTR->PFJE[fJets[i]]);
		p4s.push_back(v);
		fHemiObject.type="jet"; fHemiObject.index=fJets[i]; fHemiObject.hemi=-1;
		hemiobject.objects.push_back(fHemiObject);
	}
	for(int i=0; i< fElecs.size(); ++i){
		TLorentzVector v;
		v.SetPxPyPzE(fTR->ElPx[fElecs[i]], fTR->ElPy[fElecs[i]], fTR->ElPz[fElecs[i]],fTR->ElE[fElecs[i]]);
		p4s.push_back(v);
		fHemiObject.type="ele"; fHemiObject.index=fElecs[i]; fHemiObject.hemi=-1;
		hemiobject.objects.push_back(fHemiObject);
	}
	for(int i=0; i< fMuons.size(); ++i){
		TLorentzVector v;
		v.SetPxPyPzE(fTR->MuPx[fMuons[i]], fTR->MuPy[fMuons[i]], fTR->MuPz[fMuons[i]], fTR->MuE[fMuons[i]]);
		p4s.push_back(v);
		fHemiObject.type="muo"; fHemiObject.index=fMuons[i]; fHemiObject.hemi=-1;
		hemiobject.objects.push_back(fHemiObject);
	}
	
	// get pseudojets and minDHT
	TLorentzVector pj1, pj2;
	hemiobject.minDHT  =MinDeltaHt_pseudojets(p4s, pj1, pj2);
	
    	std::vector<double> pTs; for(unsigned i=0; i<p4s.size(); i++) pTs.push_back(p4s[i].Pt());
    	const double sumPT = accumulate( pTs.begin(), pTs.end(), double(0) );
	const TLorentzVector sumP4 = accumulate( p4s.begin(), p4s.end(), TLorentzVector() );

    	hemiobject.alphaT  = 0.5 * ( sumPT - hemiobject.minDHT ) / sqrt( sumPT*sumPT - sumP4.Perp2() );
		
	TLorentzVector pmiss;
	pmiss.SetPx(fTR->PFMETpx);
	pmiss.SetPy(fTR->PFMETpy);
	
	// fill UTM
	hemiobject.UTM          = - pmiss - pj1 - pj2;

	// fill MT2
	hemiobject.MT2   =GetMT2(pj1, 0., pj2, 0., pmiss, 0.);
	
	// fill MCT
	TVector2 pmiss_vector2;
	pmiss_vector2.Set(fTR->PFMETpx, fTR->PFMETpy);
	TLorentzVector downstream(0.,0.,0.,0.); // no particles are downstream, i.e. not selected jets are upstream. 
	hemiobject.MCT   =GetMCTcorr(pj1, pj2, downstream, pmiss_vector2);
	
	// fill pseudojets
	hemiobject.pjet1 =pj1;
	hemiobject.pjet2 =pj2;

}

// ****************************************************************************************************
double MassAnalysis::GetMCTperp(TLorentzVector p1, TLorentzVector p2, TLorentzVector P_UTM){
	double m1=p1.M();
	double m2=p2.M();
	
	TVector3 p1_T_perp  =  GetMomPerp(p1, P_UTM); 
	TVector3 p2_T_perp  =  GetMomPerp(p2, P_UTM);
	
	double E1_T_perp = sqrt( pow(m1,2) + p1_T_perp.Mag2());
	double E2_T_perp = sqrt( pow(m2,2) + p2_T_perp.Mag2());	
	
//	double MCT_perp= sqrt( pow(m1,2) + pow(m2,2) + 2 *(E1_T_perp*E2_T_perp + p1_T_perp.Dot(p2_T_perp)) );
	double MCT_perp= sqrt( pow(m1,2) + pow(m2,2) + 2 *(E1_T_perp*E2_T_perp + p1_T_perp.Mag()*p2_T_perp.Mag()*cos(p1_T_perp.Angle(p2_T_perp))) );


	// compare to Tovey Impelemntation
	TLorentzVector DTM(0.,0.,0.,0.);
	TVector2       pmiss;
	pmiss.Set(-P_UTM.Px()-p1.Px()-p2.Px(), -P_UTM.Py()-p1.Py()-p2.Py());
	double ToveyMCTperp=GetToveyMCTperp(p1, p2, DTM, pmiss);
	// if(fabs(MCT_perp - ToveyMCTperp)/ToveyMCTperp > 0.01) cout << "MCTperp" << MCT_perp << " Tovey " << ToveyMCTperp << endl;


	return MCT_perp;
}

// *****************************************************************************************************
double MassAnalysis::GetToveyMCTperp(TLorentzVector p1, TLorentzVector p2, TLorentzVector DTM, TVector2 pmiss){
	fMCT = new TMctLib;
	double MCTperp = fMCT -> mcy(p1, p2, DTM, pmiss);
	delete fMCT;
	return MCTperp;
}

// *****************************************************************************************************
double MassAnalysis::GetMCTcorr(TLorentzVector p1, TLorentzVector p2, TLorentzVector DTM, TVector2 pmiss){
	fMCT = new TMctLib;
	double MCTcorr = fMCT -> mctcorr(p1, p2, DTM, pmiss, 7000., 0.);
	delete fMCT;
	return MCTcorr;
}

// ******************************************************************************************************
double MassAnalysis::GetMCTcorr(TLorentzVector p1, TLorentzVector p2, TLorentzVector UTM){
	TVector2 pmiss;
	pmiss.Set(-UTM.Px()-p1.Px()-p2.Px(), -UTM.Py()-p1.Py()-p2.Py());
	TLorentzVector DTM(0.,0.,0.,0.);

	fMCT = new TMctLib;
	double MCTcorr = fMCT -> mctcorr(p1, p2, DTM, pmiss, 7000., 0.);
	delete fMCT;
	return MCTcorr;

}
// ******************************************************************************************************
double MassAnalysis::GetMT2perp(TLorentzVector p1, TLorentzVector p2, TLorentzVector P_UTM, double m_inv){
	TVector3 p1_T_perp =  GetMomPerp(p1, P_UTM);
	TVector3 p2_T_perp =  GetMomPerp(p2, P_UTM);
	
	// double A_T_perp   = 1/2.*(p1_T_perp.Mag()*p2_T_perp.Mag()+p1_T_perp.Dot(p2_T_perp) );
	// above implementation is problematic: A_T_perp is sometimes <0 (presicion problem) 
	// workaround: implementation below
	double A_T_perp   = 1/2.*( p1_T_perp.Mag()*p2_T_perp.Mag() * (1+cos(p1_T_perp.Angle(p2_T_perp)) ));
	double MT2_perp   = sqrt(A_T_perp) + sqrt(A_T_perp + pow(m_inv,2));
	
	return MT2_perp;
}

// *******************************************************************************************************
TVector3 MassAnalysis::GetMomPerp(TLorentzVector p, TLorentzVector P_UTM){
	TVector3 p_T(p.Px(), p.Py(), 0.);
	TVector3 pUTM_T(P_UTM.Px(), P_UTM.Py(), 0.);
	double pUTM_T_norm2=pUTM_T.Mag2();
	
	TVector3 p_T_perp = p_T - 1/pUTM_T_norm2*(p_T.Dot(pUTM_T))*pUTM_T ;	
	return p_T_perp;
}

// ****************************************************************************************************
double MassAnalysis::GetMCT(TLorentzVector p1, TLorentzVector p2){
	double m1=p1.M();
	double m2=p2.M();
	double ET1=sqrt( pow(p1.Pt(),2) +pow(m1,2) );
	double ET2=sqrt( pow(p2.Pt(),2) +pow(m2,2) );
	
	double MCTsquared=pow(m1,2) + pow(m2,2) +2*(ET1*ET2+ (p1.Px() * p2.Px()) + (p1.Py() * p2.Py()));
	
	// Tovey implementation
	fMCT = new TMctLib();
	double Tovey_MCT = fMCT -> mct(p1, p2);
	delete fMCT;
//	if(fabs(Tovey_MCT -sqrt(MCTsquared))/Tovey_MCT > 0.01) cout << "Tovey_MCT "<< Tovey_MCT << " MCT " << sqrt(MCTsquared) << endl;

	return sqrt(MCTsquared);
}

// ****************************************************************************************************
double MassAnalysis::GetAlphaT(std::vector<TLorentzVector>& p4s) {
    if(p4s.size()<2) return 0;
    
    std::vector<double> pTs; for(unsigned i=0; i<p4s.size(); i++) pTs.push_back(p4s[i].Pt());
    const std::vector<double> DHTper( DeltaSumPt_permutations(p4s) );
    
    const double mDHT = *(std::min_element( DHTper.begin(), DHTper.end() ));
    const double sumPT = accumulate( pTs.begin(), pTs.end(), double(0) );
    const TLorentzVector sumP4 = accumulate( p4s.begin(), p4s.end(), TLorentzVector() );

    return 0.5 * ( sumPT - mDHT ) / sqrt( sumPT*sumPT - sumP4.Perp2() );
}

std::vector<double> MassAnalysis::DeltaSumPt_permutations(std::vector<TLorentzVector>& p4s) {
  	std::vector<std::vector<double> > ht( 1<<(p4s.size()-1) , std::vector<double>( 2, 0.) ); // initializiert einen vector ht der size 1<<(p4s.size()-1), wobei jeder eintrag ein vector (0,0) ist
  	for(unsigned i=0; i < ht.size(); i++) {                                                  // 1<<(p4s.size()-1) = 2^(p4s.size() -1)
    	for(unsigned j=0; j < p4s.size(); j++) {
			ht [i] [(i/(1<<j))%2] += p4s[j].Pt();
    	}
  	}
  	std::vector<double> deltaHT; for(unsigned i=0; i<ht.size(); i++) deltaHT.push_back(fabs(ht[i][0]-ht[i][1]));
  	return deltaHT;
}

// ***************************************************************************************************
double MassAnalysis::MinDeltaHt_pseudojets(std::vector<TLorentzVector>& p4s, TLorentzVector& pj1, TLorentzVector& pj2){

  	std::vector<std::vector<double> > ht( 1<<(p4s.size()-1) , std::vector<double>( 2, 0.) ); // initializiert einen vector ht der size 1<<(p4s.size()-1), wobei jeder eintrag ein vector (0,0) ist
												 // 1<<(p4s.size()-1) = 2^(p4s.size() -1)
       	
  	std::vector<std::vector<double> > px( 1<<(p4s.size()-1) , std::vector<double>( 2, 0.) ); 
  	std::vector<std::vector<double> > py( 1<<(p4s.size()-1) , std::vector<double>( 2, 0.) ); 
  	std::vector<std::vector<double> > pz( 1<<(p4s.size()-1) , std::vector<double>( 2, 0.) ); 
  	std::vector<std::vector<double> > E ( 1<<(p4s.size()-1) , std::vector<double>( 2, 0.) ); 

  	for(unsigned i=0; i < ht.size(); i++) {                                             
    	for(unsigned j=0; j < p4s.size(); j++) {
			ht [i] [(i/(1<<j))%2] += p4s[j].Pt();
			px [i] [(i/(1<<j))%2] += p4s[j].Px();
			py [i] [(i/(1<<j))%2] += p4s[j].Py();
			pz [i] [(i/(1<<j))%2] += p4s[j].Pz();
			E  [i] [(i/(1<<j))%2] += p4s[j].E();
    	}
  	}
  	std::vector<double> deltaHT; for(unsigned i=0; i<ht.size(); i++) deltaHT.push_back(fabs(ht[i][0]-ht[i][1]));
	const double mDHT = *(std::min_element( deltaHT.begin(), deltaHT.end() ));
	int pos=distance(deltaHT.begin(), min_element(deltaHT.begin(), deltaHT.end()));
	pj1.SetPxPyPzE(px[pos][0], py[pos][0], pz[pos][0], E[pos][0]);
	pj2.SetPxPyPzE(px[pos][1], py[pos][1], pz[pos][1], E[pos][1]);
	return mDHT;
}


// ****************************************************************************************************
double MassAnalysis::GetMT2(TLorentzVector v1, double mv1, TLorentzVector v2, double mv2, TLorentzVector p_unobs, int m_invisible){
	double pa[3];
	double pb[3];
	double pmiss[3];

	pmiss[1]=p_unobs.Px();
	pmiss[2]=p_unobs.Py();

	pa[0]=mv1;
	pa[1]=v1.Px();
	pa[2]=v1.Py();

	pb[0]=mv2;
	pb[1]=v2.Px();
	pb[2]=v2.Py();
	
	fMT2 = new Davismt2();
	fMT2->set_momenta(pa, pb, pmiss);
	fMT2->set_mn(m_invisible);
	double MT2=fMT2->get_mt2();
	delete fMT2;
	return MT2;
}


// ****************************************************************************************************
vector<TLorentzVector> MassAnalysis::GetLepton4Momenta(){
	vector<TLorentzVector> momenta;
	for(int i=0; i<fElecs.size(); ++i){
		TLorentzVector p;
		p.SetPtEtaPhiE(fTR->ElPt[fElecs[i]],fTR->ElEta[fElecs[i]], fTR->ElPhi[fElecs[i]], fTR->ElE[fElecs[i]]);
		momenta.push_back(p);
	}
	for(int i=0; i<fMuons.size(); ++i){
		TLorentzVector p;
		p.SetPtEtaPhiM(fTR->MuPt[fMuons[i]],fTR->MuEta[fMuons[i]], fTR->MuPhi[fMuons[i]],0.105);
		momenta.push_back(p);
	}
	return momenta;
}


// ****************************************************************************************************
void MassAnalysis::End(){
	// cout interesting events
	cout << " *************************************************************** " << endl;
	cout << " interesting events                                              " << endl;
	for(int i=0; i< interesting_Event.size(); ++i){
		cout << interesting_Type[i] << "Run:Lumi:Evt" << " " << interesting_Run[i] << ":" << interesting_Lumi[i]
		     << ":" << interesting_Event[i] << " value " << interesting_value[i] << endl;	     
	}
	cout << " *************************************************************** " << endl;
	
	fHistFile->cd();	

	// write tree
	fATree->Write();
	fHistFile                ->Close();
}
