#include "TreeCleaner.hh"
#include "helper/Utilities.hh"

using namespace std;

TreeCleaner::TreeCleaner(TreeReader *tr) : UserAnalysisBase(tr){
	SetClean(true);
	SetSkim(true);
}

TreeCleaner::~TreeCleaner(){
}

void TreeCleaner::Begin(){
	ReadCleaningParameters();
	StatInit();
	Reset();
	if(fSkim){
		fCleanTreeFile = new TFile(fOutputDir + "CleanTree.root", "RECREATE");
		fCleanTreeFile->mkdir("analyze", "analyze");
		fCleanTreeFile->cd("analyze");
		fCleanTree = fTR->fChain->CloneTree(0);
		fCleanTree->CopyAddresses(fTR->fChain);
	}
}

void TreeCleaner::Analyze(){
// Performs the tagging, isolation and the object cleaning

	DoTagging();
	if(fClean) DoCleaning();
	if(fSkim)  DoSkimTree();
}

void TreeCleaner::DoTagging(){
// Performs only the tagging and isolation
	InitCleaning();
	TagCleanObjects();
	DecideIso();
	TagDuplObjects();
}

void TreeCleaner::DoCleaning(){
// Performs only the object cleaning
//   assumes that DoTagging() has been called
	DoCleanObjects();
	StatFill();
}

void TreeCleaner::DoSkimTree(){
// Performs only the tree skimming, removing bad objects
//   assumes that DoTagging() and DoCleaning() have been called

	if( fClean ){
		// Clean and skim the tree
		int iclean = 0;
		int ntkmu(0), ngbmu(0);
		for( int i = 0; i < fTR->NMus; i++ ){
			if( fTR->MuGood[i] != 0 ) continue;
			PutMuon(iclean, i);
			iclean++;
			if(fTR->MuIsTrackerMuon[i]) ntkmu++;
			if(fTR->MuIsGlobalMuon[i]) ngbmu++;
		}
		fTR->NMus = iclean;
		fTR->NTMus = ntkmu;
		fTR->NGMus = ngbmu;

		iclean = 0;
		for( int i = 0; i < fTR->NEles; i++ ){
			if( fTR->ElGood[i] != 0 ) continue;
			PutElectron(iclean, i);
			iclean++;
		}
		fTR->NEles = iclean;

		iclean = 0;
		for( int i = 0; i < fTR->NPhotons; i++ ){
			if( fTR->PhoGood[i] != 0 ) continue;
			PutPhoton(iclean, i);
			iclean++;
		}
		fTR->NPhotons = iclean;

		iclean = 0;
		for( int i = 0; i < fTR->NJets; i++ ){
			if( fTR->JGood[i] != 0 ) continue;
			PutJet(iclean, i);
			iclean++;
		}
		fTR->NJets = iclean;

		// At this point the arrays only contain clean objects
		// The non-clean object are overwritten
		fCleanTree->Fill();
	}
}

void TreeCleaner::End(){
	if(fClean){
		if(fVerbose) StatPrint();
		StatHistos();
	}

	if(fSkim){		
		fCleanTreeFile->cd("analyze");
		fCleanTree->Write();
		fCleanTreeFile->Close();
	}
	fHstatFile->cd();
	fHstatHistos->Write();
	fHstatFile->Close();
	
}

void TreeCleaner::Reset(){
	fR12 = -999;
	fR21 = -999;
}

void TreeCleaner::ReadCleaningParameters(const char* filename){
	ifstream IN(filename);
	char buffer[200];
	char ParName[100];
	float ParValue;

	bool verbose(true);
	bool ok(false);

	while( IN.getline(buffer, 200, '\n') ){
		ok = false;
		if (buffer[0] == '#') {continue;} // Skip lines commented with '#'
		sscanf(buffer, "%s %f", ParName, &ParValue);

		// final acceptance cuts
		if( !strcmp(ParName, "minJetPt") ){
			fMinJetPt = float(ParValue); ok = true;
		}

		// -- Primary vertex:
		if( !strcmp(ParName, "chisqVxmax") ){
			fClean_chisqVxmax = float(ParValue); ok = true;
		}
		if( !strcmp(ParName, "dRVxmax") ){
			fClean_dRVxmax = float(ParValue); ok = true;
		}
		if( !strcmp(ParName, "dzVxmax") ){
			fClean_dzVxmax = float(ParValue); ok = true;
		}
		if( !strcmp(ParName, "sumPtTkfromVxmin") ){
			fClean_sumPtTkfromVxmin = float(ParValue); ok = true;
		}
		if( !strcmp(ParName, "PrimVtxNdofmin") ){
			fClean_PrimVtxNdofmin = int(ParValue); ok = true;
		}
		
		if( !strcmp(ParName, "distVxmax") ){
			fClean_distVxmax = float(ParValue); ok = true;
		}
		

		// -- Muons:
		if( !strcmp(ParName, "MuonDPbyPmax") ){
			fClean_MuonDPbyPmax = float(ParValue); ok = true;
		}
		if( !strcmp(ParName, "MuonChi2max") ){
			fClean_MuonChi2max = float(ParValue); ok = true;
		}
		if( !strcmp(ParName, "MuonNHitsmin") ){
			fClean_MuonNHitsmin = float(ParValue); ok = true;
		}
		if( !strcmp(ParName, "dRSSmuonmax") ){
			fClean_dRSSmuonmax = float(ParValue); ok = true;
		}

		// -- Electrons:
		if( !strcmp(ParName, "ElecHoverEBarmax") ){
			fClean_ElecHoverEBarmax = float(ParValue); ok = true;
		}
		if( !strcmp(ParName, "ElecHoverEEndmax") ){
			fClean_ElecHoverEEndmax = float(ParValue); ok = true;
		}
		if( !strcmp(ParName, "ElecSigmaEtaEtaBarmin") ){
			fClean_ElecSigmaEtaEtaBarmin = float(ParValue); ok = true;
		}
		if( !strcmp(ParName, "ElecSigmaEtaEtaBarmax") ){
			fClean_ElecSigmaEtaEtaBarmax = float(ParValue); ok = true;
		}
		if( !strcmp(ParName, "ElecSigmaEtaEtaEndmax") ){
			fClean_ElecSigmaEtaEtaEndmax = float(ParValue); ok = true;
		}
		if( !strcmp(ParName, "ElecEoverPInBarmin") ){
			fClean_ElecEoverPInBarmin = float(ParValue); ok = true;
		}
		if( !strcmp(ParName, "ElecEoverPInEndmin") ){
			fClean_ElecEoverPInEndmin = float(ParValue); ok = true;
		}
		if( !strcmp(ParName, "ElecDeltaEtaInBarmax") ){
			fClean_ElecDeltaEtaInBarmax = float(ParValue); ok = true;
		}
		if( !strcmp(ParName, "ElecDeltaEtaInEndmax") ){
			fClean_ElecDeltaEtaInEndmax = float(ParValue); ok = true;
		}
		if( !strcmp(ParName, "ElecDeltaPhiInBarmax") ){
			fClean_ElecDeltaPhiInBarmax = float(ParValue); ok = true;
		}
		if( !strcmp(ParName, "ElecDeltaPhiInEndmax") ){
			fClean_ElecDeltaPhiInEndmax = float(ParValue); ok = true;
		}
		if( !strcmp(ParName, "ElecDeltaPhiOutBarmax") ){
			fClean_ElecDeltaPhiOutBarmax = float(ParValue); ok = true;
		}
		if( !strcmp(ParName, "ElecDeltaPhiOutEndmax") ){
			fClean_ElecDeltaPhiOutEndmax = float(ParValue); ok = true;
		}
		if( !strcmp(ParName, "ElecConvPartTrackDistmax") ){
			fClean_ElecConvPartTrackDistmax = float(ParValue); ok = true;
		}
		if( !strcmp(ParName, "ElecConvPartTrackDCotmax") ){
			fClean_ElecConvPartTrackDCotmax = float(ParValue); ok = true;
		}
		if( !strcmp(ParName, "ElecNMissHitsmax") ){
			fClean_ElecNMissHitsmax = float(ParValue); ok = true;
		}
		if( !strcmp(ParName, "dRSSelecmax") ){
			fClean_dRSSelecmax = float(ParValue); ok = true;
		}

		// -- Photons:
		if( !strcmp(ParName, "PhotHoverEBarmax") ){
			fClean_PhotHoverEBarmax = float(ParValue); ok = true;
		}
		if( !strcmp(ParName, "PhotHoverEEndmax") ){
			fClean_PhotHoverEEndmax = float(ParValue); ok = true;
		}
		if( !strcmp(ParName, "PhotSigmaEtaEtaBarmax") ){
			fClean_PhotSigmaEtaEtaBarmax = float(ParValue); ok = true;
		}
		if( !strcmp(ParName, "PhotSigmaEtaEtaEndmax") ){
			fClean_PhotSigmaEtaEtaEndmax = float(ParValue); ok = true;
		}
		if( !strcmp(ParName, "PhotSigmaEtaEtaBarmin") ){
			fClean_PhotSigmaEtaEtaBarmin = float(ParValue); ok = true;
		}

		// -- Jets:
		if( !strcmp(ParName, "FracEmmaxJet") ){
			fClean_FracEmmaxJet = float(ParValue); ok = true;
		}
		if( !strcmp(ParName, "FracEmminJet") ){
			fClean_FracEmminJet = float(ParValue); ok = true;
		}
		if( !strcmp(ParName, "FracChminJet") ){
			fClean_FracChminJet = float(ParValue); ok = true;
		}
		if( !strcmp(ParName, "JID_n90Hitsmin") ){
			fClean_JID_n90Hitsmin = float(ParValue); ok = true;
		}
		if( !strcmp(ParName, "JID_HPDmax") ){
			fClean_JID_HPDmax = float(ParValue); ok = true;
		}
		if( !strcmp(ParName, "JID_RBXmax") ){
			fClean_JID_RBXmax = float(ParValue); ok = true;
		}
		if( !strcmp(ParName, "deltaRElecJetmax") ){
			fClean_deltaRElecJetmax = float(ParValue); ok = true;
		}
		if( !strcmp(ParName, "elecbyJetEratio") ){
			fClean_elecbyJetEratio = float(ParValue); ok = true;
		}

		// -- Isolation:
		if( !strcmp(ParName, "MuonIsomax") ){
			fClean_MuonIsomax = float(ParValue); ok = true;
		}
		if( !strcmp(ParName, "ElecIsomax") ){
			fClean_ElecIsomax = float(ParValue); ok = true;
		}
		if( !strcmp(ParName, "PhotIsomax") ){
			fClean_PhotIsomax = float(ParValue); ok = true;
		}

		// -- Event Cleaning:
		if( !strcmp(ParName, "FracChmin") ){
			fClean_FracChmin = float(ParValue); ok = true;
		}
		if( !strcmp(ParName, "FracEmmin") ){
			fClean_FracEmmin = float(ParValue); ok = true;
		}
		if( !strcmp(ParName, "JetBadHardPtmin") ){
			fClean_JetBadHardPtmin = float(ParValue); ok = true;
		}

		// -- MET:
		if( !strcmp(ParName, "METmin") ){
			fClean_METmin = float(ParValue); ok = true;
		}
		if( !strcmp(ParName, "dPhiJetMETmin") ){
			fClean_dPhiJetMETmin = float(ParValue); ok = true;
		}
		if( !strcmp(ParName, "dR12min") ){
			fClean_dR12min = float(ParValue); ok = true;
		}
		if( !strcmp(ParName, "dR21min") ){
			fClean_dR21min = float(ParValue); ok = true;
		}

		if(!ok) cout << "%% TreeCleaner::ReadCleaningParameters ==> ERROR: Unknown variable " << ParName << endl;
	}
}

void TreeCleaner::TagCleanObjects(void){
// Steering for tagging clean/bad objects

	// Primary vertex
	fTR->PrimVtxGood = CleanPrimaryVertex();

	// Muons
	for( int ichk = 0; ichk < fTR->NMus; ++ichk ){
		fTR->MuGood[ichk] = 10*IsFromPrimaryVx(1, ichk);
		fTR->MuGood[ichk] += CleanMuon(ichk);
	}

	// Electrons
	for( int ichk = 0; ichk < fTR->NEles; ++ichk ){
		fTR->ElGood[ichk] = 10*IsFromPrimaryVx(2, ichk);
		fTR->ElGood[ichk] += CleanElectron(ichk);
	}

	// Photons
	for( int ichk = 0; ichk < fTR->NPhotons; ++ichk ){
		fTR->PhoGood[ichk] = CleanPhoton(ichk);
	}

	// Jets
	for( int ichk = 0; ichk < fTR->NJets; ++ichk ){
		fTR->JGood[ichk] = 10*IsFromPrimaryVx(4, ichk);
		fTR->JGood[ichk] += CleanJet(ichk, 1);
	}
	return;
}

void TreeCleaner::TagDuplObjects(void){
// Steering for tagging clean/bad objects

	// Duplication (only after cleanness has been checked)
	int jbad = 0;
	for( int ichk = 0; ichk < fTR->NMus; ++ichk )  if( DuplicateMuon(ichk)     ) fTR->MuGood[ichk] += 100;
	for( int ichk = 0; ichk < fTR->NEles; ++ichk ) if( DuplicateElectron(ichk) ) fTR->ElGood[ichk] += 100;
	for (int ichk = 0; ichk < fTR->NPhotons; ++ichk) if ( DuplPhotonElectron(ichk) ) fTR->PhoGood[ichk] += 100;
	//	for( int ichk = 0; ichk < fTR->NJets; ++ichk ) {
	//  should be done later
	//		if( ElectronJet(ichk)     ) fTR->JGood[ichk]  += 100;
	//		if( PhotonJet(ichk)       ) fTR->JGood[ichk]  += 200;
	//		if (fTR->JGood[ichk]%10 > 0 && fTR->JPt[ichk] > fClean_JetBadHardPtmin) jbad = 4;
	//	}

	// Event and MET can only be done later (after bad objects are removed)
	fTR->GoodEvent = jbad;
	return;
}

int TreeCleaner::CleanPrimaryVertex(void){
// Verifies the primary vertex quality
// returns iBad = 1 for no charged tracks
//              = 2 for bad normalized chi squared AND 1
//              = 3 for incompatible with beamspot AND 1 AND 2
//              = 4 for insufficient track pT AND 1 AND 2 AND 3
	int iBad = 0;
// Check that there are tracks at the Primary Vertex
	if( fTR->PrimVtxNdof < fClean_PrimVtxNdofmin ){
		iBad = 1;
		return iBad;
	}

// Check the chi2/ndof
	if( fTR->PrimVtxNChi2 > fClean_chisqVxmax || fTR->PrimVtxNChi2 < 0. ){
		iBad = 2;
		return iBad;
	}

// Check compatibility of vertex with beam spot
	double xVx = fTR->PrimVtxx - fTR->Beamspotx;
	double yVx = fTR->PrimVtxy - fTR->Beamspoty;
	double zVx = fTR->PrimVtxz - fTR->Beamspotz;
	double rVx = sqrt(xVx*xVx + yVx*yVx);
	if (rVx > fClean_dRVxmax || fabs(zVx) > fClean_dzVxmax) {
		iBad = 3;
		return iBad;
	}

// Check that there is sufficient Et in the tracks
	if( fTR->PrimVtxPtSum < fClean_sumPtTkfromVxmin ){
		iBad = 4;
		return iBad;
	}
	return iBad;
}

int TreeCleaner::IsFromPrimaryVx(int ipart, int ichk){
// Checks whether the object is compatible with the primary vertex
// ipart = 1 for muon
//       = 2 for electron
//       = 3 for photon
//       = 4 for jet
// returns iBad = 1 for incompatible with primary vertex
	int iBad = 0;
	
// take error from vertex and from track extrapolation into account
	double drVxsq = fTR->PrimVtxxE*fTR->PrimVtxxE + fTR->PrimVtxyE*fTR->PrimVtxyE;
	double d0, dd0, dz, ddz;
	if( ipart <= 0 || ichk < 0 ){
		return iBad;
	} else if (fTR->PrimVtxGood != 0) {
		return iBad;
	} else if( ipart == 1 ){ // Muons
		d0  = fTR->MuD0PV[ichk];
		dd0 = sqrt(fTR->MuD0E[ichk]*fTR->MuD0E[ichk] + drVxsq);
		dz  = fTR->MuDzPV[ichk];
		ddz = sqrt(fTR->MuDzE[ichk]*fTR->MuDzE[ichk] + fTR->PrimVtxzE*fTR->PrimVtxzE);
	} else if( ipart == 2 ){ // Electrons
		d0  = fTR->ElD0PV[ichk];
		dd0 = sqrt(fTR->ElD0E[ichk]*fTR->ElD0E[ichk] + drVxsq);
		dz  = fTR->ElDzPV[ichk];
		ddz = sqrt(fTR->ElDzE[ichk]*fTR->ElDzE[ichk] + fTR->PrimVtxzE*fTR->PrimVtxzE);
	} else if( ipart == 3 ){ // Photons
		return 1; // Photons not implemented yet
	} else if( ipart == 4 ){ // Jets
		d0  = 0.;
		dd0 = 0.001;
		dz  = 0.;
		ddz = 0.001;
//		if (fTR->JVtxNChi2 > 0) { // not defined when bad
		double dztry = fTR->JVtxz[ichk] - fTR->PrimVtxz;
		if (dztry > -100.) {
			d0  = sqrt((fTR->JVtxx[ichk]-fTR->PrimVtxx)*(fTR->JVtxx[ichk]-fTR->PrimVtxx)
				+ (fTR->JVtxy[ichk]-fTR->PrimVtxy)*(fTR->JVtxy[ichk]-fTR->PrimVtxy) );
			dd0 = sqrt(fTR->JVtxExx[ichk] + fTR->JVtxEyy[ichk] + drVxsq);
			dz  = fTR->JVtxz[ichk] - fTR->PrimVtxz;
			ddz = sqrt(fTR->JVtxEzz[ichk] + fTR->PrimVtxzE*fTR->PrimVtxzE);
		}
	}

// test that the distance is not too large
	if( fabs(d0) > fClean_distVxmax * dd0 || fabs(dz) > fClean_distVxmax * ddz ){
		iBad = 1;
	}
	return iBad;
}

int TreeCleaner::CleanMuon(int ichk){
// Verifies the muon identification quality
// returns iBad = 0 for good muons
//              = 1 for bad Delta pT / pT
//              = 2 for bad normalized chi squared
//              = 3 for too few valid hits in tracker
	int iBad = 0;
	if( ichk < 0 ) return -1;
	if( fTR->MuPtE[ichk] >= fClean_MuonDPbyPmax * fTR->MuPt[ichk] )  iBad = 1; // Maximum Delta p / p
	else if( fTR->MuNChi2[ichk] > fClean_MuonChi2max )               iBad = 2; // Maximum Chisquared
	else if( fTR->MuNTkHits[ichk] < fClean_MuonNHitsmin )            iBad = 3; // Minimum number of valid hits
	if ( !fTR->MuIsGlobalMuon[ichk] || !fTR->MuIsTrackerMuon[ichk])  iBad = 4;
	return iBad;
}

bool TreeCleaner::DuplicateMuon(int ichk){
// Checks for duplicate muons
	// fClean_dRSSmuonmax = 0.1;

	if( ichk < 0 ) return false;

	for( int j = 0; j < fTR->NMus; ++j ){
		if( j == ichk ) continue;
		if( fTR->MuCharge[ichk] != fTR->MuCharge[j] ) continue;

		double deltaR = Util::GetDeltaR(fTR->MuEta[ichk], fTR->MuEta[j], fTR->MuPhi[ichk], fTR->MuPhi[j]);
		if( deltaR > fClean_dRSSmuonmax ) continue;

		// Both are bad or both are good -> compare them
		if( (fTR->MuGood[ichk] == 0 && fTR->MuGood[j] == 0) || (fTR->MuGood[ichk] != 0 && fTR->MuGood[j] != 0) ){
			if( fTR->MuPtE[ichk]/fTR->MuPt[ichk] >= fTR->MuPtE[j]/fTR->MuPt[j] ) return true;
		}
		// One is good, one is bad -> take the good one
		else if( fTR->MuGood[j] == 0 && fTR->MuGood[ichk] != 0 ) return true;
	}
	return false;
}

int TreeCleaner::CleanElectron(int ichk){
// Verifies the electron identification quality
// returns iBad = 0 for good muons
//              = 1 for bad H/E
//              = 2 for bad shower shape
//              = 3 for bad matching of Ecal and track
	bool useHoverE      = true;
	bool useSigmaEtaEta = true;
	bool useEoverPIn    = true;
	bool useDeltaEtaIn  = true;
	bool useDeltaPhiIn  = true;
	bool useDeltaPhiOut = true;
	if( ichk < 0 ) return -1;

	double hOverE      = fTR->ElHcalOverEcal[ichk];
	double sigmaee     = fTR->ElSigmaIetaIeta[ichk];
	double eOverPin    = fTR->ElESuperClusterOverP[ichk];
	double deltaEtaIn  = fTR->ElDeltaEtaSuperClusterAtVtx[ichk];
	double deltaPhiIn  = fTR->ElDeltaPhiSuperClusterAtVtx[ichk];
	double deltaPhiOut = fTR->ElDeltaPhiSeedClusterAtCalo[ichk];
	double nMissHits   = fTR->ElNumberOfMissingInnerHits[ichk];
	double trackDist   = fTR->ElConvPartnerTrkDist[ichk];
	double trackDCot   = fTR->ElConvPartnerTrkDCot[ichk]; 

	// Distinguish between barrel and endcap electrons
	if( fabs(fTR->ElEta[ichk]) < 1.479 ){ // Barrel
		if(useHoverE)      if( hOverE            > fClean_ElecHoverEBarmax      ) return 1;
		if(useSigmaEtaEta) if( sigmaee           > fClean_ElecSigmaEtaEtaBarmax ) return 2;
		if(useEoverPIn)    if( eOverPin          < fClean_ElecEoverPInBarmin    ) return 3;
		if(useDeltaEtaIn)  if( fabs(deltaEtaIn)  > fClean_ElecDeltaEtaInBarmax  ) return 3;
		if(useDeltaPhiIn)  if( fabs(deltaPhiIn)  > fClean_ElecDeltaPhiInBarmax  ) return 3;
		if(useDeltaPhiOut) if( fabs(deltaPhiOut) > fClean_ElecDeltaPhiOutBarmax ) return 3;
		if(useSigmaEtaEta) if( sigmaee           < fClean_ElecSigmaEtaEtaBarmin ) return 4;
	}
	else{ // EndCap
		if(useHoverE)      if( hOverE            > fClean_ElecHoverEEndmax      ) return 1;
		if(useSigmaEtaEta) if( sigmaee           > fClean_ElecSigmaEtaEtaEndmax ) return 2;
		if(useEoverPIn)    if( eOverPin          < fClean_ElecEoverPInEndmin    ) return 3;
		if(useDeltaEtaIn)  if( fabs(deltaEtaIn)  > fClean_ElecDeltaEtaInEndmax  ) return 3;
		if(useDeltaPhiIn)  if( fabs(deltaPhiIn)  > fClean_ElecDeltaPhiInEndmax  ) return 3;
		if(useDeltaPhiOut) if( fabs(deltaPhiOut) > fClean_ElecDeltaPhiOutEndmax ) return 3;
	}
	if (  nMissHits > fClean_ElecNMissHitsmax )         return 5;
	if (  trackDist < fClean_ElecConvPartTrackDistmax && trackDCot < fClean_ElecConvPartTrackDCotmax ) return 6;
	return 0;
}

bool TreeCleaner::DuplicateElectron(int ichk){
// Checks for duplicate electrons
	if( ichk < 0 ) return false;

	int j = fTR->ElDuplicateEl[ichk];
	if( j < 0 ) return false;

	if( Util::GetDeltaR(fTR->ElEta[ichk], fTR->ElEta[j], fTR->ElPhi[ichk], fTR->ElPhi[j]) < fClean_dRSSelecmax ){
		// Both are bad or both are good -> compare them
		if( (fTR->ElGood[ichk] == 0 && fTR->ElGood[j] == 0) || (fTR->ElGood[ichk] != 0 && fTR->ElGood[j] != 0) ){
			double elecEoP = fTR->ElESuperClusterOverP[ichk];
			double newEoP  = fTR->ElESuperClusterOverP[j];
			if( fabs(fTR->ElESuperClusterOverP[ichk]-1.) > fabs(fTR->ElESuperClusterOverP[j]-1.) ) return true;
		}
		// One is good, one is bad -> take good one
		else if( fTR->ElGood[j] == 0 && fTR->ElGood[ichk] != 0 ) return true;
	}
	return false;
}

int TreeCleaner::CleanPhoton(int ichk){
// Verifies the photon identification quality
// returns iBad = 1 for bad H/E
// Still to be completed ****
	if( ichk < 0 ) return 0;
	bool useHoverE      = true;
	bool useSigmaEtaEta = true;

	// Distinguish between barrel and endcap photons
	if( fabs(fTR->PhoEta[ichk]) < 1.479 ){ // Barrel
		if(useHoverE)      if( fTR->PhoHoverE[ichk]        > fClean_PhotHoverEBarmax      ) return 1;
		if(useSigmaEtaEta) if( fTR->PhoSigmaIetaIeta[ichk] > fClean_PhotSigmaEtaEtaBarmax ) return 2;
		if(useSigmaEtaEta) if( fTR->PhoSigmaIetaIeta[ichk] < fClean_PhotSigmaEtaEtaBarmin ) return 3;
	}
	else{ // EndCap
		if(useHoverE)      if( fTR->PhoHoverE[ichk]        > fClean_PhotHoverEEndmax      ) return 1;
		if(useSigmaEtaEta) if( fTR->PhoSigmaIetaIeta[ichk] > fClean_PhotSigmaEtaEtaEndmax ) return 2;
	}

	return 0;
}

bool TreeCleaner::DuplPhotonElectron(int ichk){
// Checks for duplicate photons with electrons
	if( ichk < 0 ) return false;

	int j = fTR->PhoIsElDupl[ichk];
/*	// fiddle because duplicate flag not in ntuple
	double DR = 999.;
	int j = -1;
	for (int i = 0; i < fTR->NEles; ++i) {
	    double DRnew =  Util::GetDeltaR(fTR->PhoEta[ichk], fTR->ElEta[i], fTR->PhoPhi[ichk], fTR->ElPhi[i]);
	    if (DRnew < DR) {
	      DR = DRnew;
	      j = i;
	    }
	}
	if( j < 0 || DR > 0.3) return false;
*/
	// Both are bad or both are good -> photon is duplicate
	if( (fTR->PhoGood[ichk] == 0 && fTR->ElGood[j] == 0) || (fTR->PhoGood[ichk] != 0 && fTR->ElGood[j] != 0) ) return true;

	// One is good, one is bad -> take good one
      	else if( fTR->ElGood[j] == 0 && fTR->PhoGood[ichk] != 0 ) return true;
	return false;
}

int TreeCleaner::CleanJet(int ichk, int itype){
// Verifies the jet reconstruction quality
// and flags jets made from electrons
// (electrons should be filled before cleaning the jets)
// returns iBad = 0 for good jets
//              = 1 Et < Pt
//              = 2 for too large EM fraction
//              = 3 for too small EM fraction
//              = 4 for bad Trk pT fraction
//              = 5 for too small JID_n90Hits
//              = 6 for too high JID_HPD (= fraction in single HPD)
//              = 7 for too high JID_RBX
	int iBad = 0;
	if( ichk < 0 ) return 0;

// veto jets with E<p
	if( fTR->JEt[ichk] - fTR->JPt[ichk] < -0.0001 ) return 1;

// check for noise
	if (itype == 1 || itype == 3) {
	    if( fTR->JID_n90Hits[ichk] < fClean_JID_n90Hitsmin ) return 5;
	    if( fTR->JID_HPD[ichk] > fClean_JID_HPDmax ) return 6;
	    if (fTR->JID_RBX[ichk] > fClean_JID_RBXmax ) return 7;
	}

// check EM and track pT fractions
	// the fEM and fCh checks are better done afterwards
	if (itype == 2 || itype == 3) {
	    if( fTR->JEMfrac[ichk] > fClean_FracEmmaxJet ) return 2;
	    if( fTR->JEMfrac[ichk] < fClean_FracEmminJet ) return 3;
	    double fChmin = fClean_FracChminJet;
	    if (fabs(fTR->JEta[ichk]) > 1.9) fChmin = fClean_FracChminJet * (1. - fabs(fTR->JEta[ichk]) + 1.9);
	    fChmin < 0. ? 0. : fChmin;
	    if( fTR->JChfrac[ichk] < fChmin && fabs(fTR->JEta[ichk]) < 2.9) return 4;
	}
	return 0;
}

bool TreeCleaner::ElectronJet(int ichk){
// checks for jets made from electrons
// ichk = index of the jet
	if( ichk < 0 ) return false;

	bool isDuplicate = false;

// veto jets made of electrons
	for( int j = 0; j < fTR->NEles; ++j ){
		if( fTR->ElIsInJet[j] < 0 ) continue;
		if( fTR->ElIsInJet[j] != ichk ) continue;
		if (fTR->ElGood[j] != 0 || fTR->ElIsIso[j] == 0) continue;
		if( Util::GetDeltaR(fTR->JEta[ichk], fTR->ElEta[j], fTR->JPhi[ichk], fTR->ElPhi[j]) > fClean_deltaRElecJetmax ) continue;
		double jetEuncorr = fTR->JE[ichk] / fTR->JEcorr[ichk];
		if( fTR->ElSharedEnergy[j] > fClean_elecbyJetEratio * jetEuncorr ){
			isDuplicate = true;
			break;
		}
	}
	return isDuplicate;
}

bool TreeCleaner::PhotonJet(int ichk){
// checks for jets made from photons
// ichk = index of the jet
	if( ichk < 0 ) return false;

	bool isDuplicate = false;

// veto jets made of photons
	for( int j = 0; j < fTR->NPhotons; ++j ){
		if( fTR->PhoIsInJet[j] < 0 ) continue;
		if( fTR->PhoIsInJet[j] != ichk ) continue;
		if (fTR->PhoGood[j] != 0 || fTR->PhoIsIso[j] == 0) continue;
		if( Util::GetDeltaR(fTR->JEta[ichk], fTR->PhoEta[j], fTR->JPhi[ichk], fTR->PhoPhi[j]) > fClean_deltaRElecJetmax ) continue;
		double jetEuncorr = fTR->JE[ichk] / fTR->JEcorr[ichk];
		if( fTR->PhoSharedEnergy[j] > fClean_elecbyJetEratio * jetEuncorr ){
			isDuplicate = true;
			break;
		}
	}
	return isDuplicate;
}

void TreeCleaner::DecideIso(void){
// decide whether objects are isolated or not
	// Muons
	for( int ichk = 0; ichk < fTR->NMus; ++ichk ){
		if( fTR->MuRelIso03[ichk] < fClean_MuonIsomax ) fTR->MuIsIso[ichk] = 1;
		else fTR->MuIsIso[ichk] = 0;
	}

	// Electrons
	for( int ichk = 0; ichk < fTR->NEles; ++ichk ){
		if( fTR->ElRelIso04[ichk] < fClean_ElecIsomax ) fTR->ElIsIso[ichk] = 1;
		else fTR->ElIsIso[ichk] = 0;
	}

	// Photons
	for ( int ichk = 0; ichk < fTR->NPhotons; ++ichk){
		if ( fTR->PhoIso03[ichk] < fClean_PhotIsomax ) fTR->PhoIsIso[ichk] = 1;
		else fTR->PhoIsIso[ichk] = 0;
	}

	return;
}

void TreeCleaner::InitCleaning(){
// Initializes the object cleaning
// Has to be called once each event!
	if(fClean){
		fNMuClean = 0;
		fNElClean = 0;
		fNJClean  = 0;
	} else {
		fNMuClean = fTR->NMus;
		fNElClean = fTR->NEles;
		fNJClean  = fTR->NJets;
	}
	return;
}

void TreeCleaner::DoCleanObjects(void){
// Steering for removal of bad objects
// Should be initialized by calling InitCleaning
// Should only be called after deciding on isolation

	// Check the primary vertex
	if( fTR->PrimVtxGood != 0 ){
		fNMuClean = 0;
		fNElClean = 0;
		fNJClean = 0;
	}

	// Muons
	fNMuClean = 0;
	for( int ichk = 0; ichk < fTR->NMus; ++ichk ){
		if( fTR->MuGood[ichk] != 0 || fTR->MuIsIso[ichk] != 1 ) continue;
		fNMuClean++;
	}

	// Electrons
	fNElClean = 0;
	for( int ichk = 0; ichk < fTR->NEles; ++ichk ){
		if( fTR->ElGood[ichk] != 0 || fTR->ElIsIso[ichk] != 1 ) continue;
		fNElClean++;
	}

	// Photons
	fNPhClean = 0;
	for( int ichk = 0; ichk < fTR->NPhotons; ++ichk ){
		if( fTR->PhoGood[ichk] != 0 || fTR->PhoIsIso[ichk] != 1 ) continue;
		fNPhClean++;
	}

	// Clean the jets
	fNJClean = 0;
	for( int ichk = 0; ichk < fTR->NJets; ++ichk ){
		if( fTR->JGood[ichk] != 0 ) continue;
		fNJClean++;
		for( int imu = 0; imu < fTR->NMus; ++imu ){
		        if( fTR->MuGood[imu] != 0 ) continue;
		        if ( ! fTR->MuIsIso[imu] ) {
			         int ij = FindNearestJet(fTR->MuEta[imu], fTR->MuPhi[imu]);
				 double DR = Util::GetDeltaR(fTR->MuEta[imu], fTR->JEta[ichk], fTR->MuPhi[imu], fTR->JPhi[ichk]);
				 if (DR > fClean_deltaRElecJetmax) continue;
				 if (ij != ichk) continue;
			         AddToJet(1, imu, ichk);
		                 fTR->MuGood[imu] += 1000;
			}
		}
		for( int iel = 0; iel < fTR->NEles; ++iel ){
			if( fTR->ElIsInJet[iel] != ichk ) continue;
			// ignore duplicated electrons
      			if (fTR->ElGood[iel] > 99) continue;
			if( fTR->ElIsIso[iel] ) {
			        SubtrFromJet(2, iel, ichk);
			        if (fTR->JPt[ichk] < fMinJetPt) fTR->JGood[ichk] += 100;
			}
			else {
			        AddToJet(2, iel, ichk);
				fTR->ElGood[iel] += 1000;
			}
		}
		for( int iph = 0; iph < fTR->NPhotons; ++iph ){
			if( fTR->PhoIsInJet[iph] != ichk ) continue;
			// ignore duplicated photons
      			if (fTR->PhoGood[iph] > 99) continue;
			if( fTR->PhoIsIso[iph] ) {
			        SubtrFromJet(3, iph, ichk);
			        if (fTR->JPt[ichk] < fMinJetPt) fTR->JGood[ichk] += 200;
			}
			else {
			        AddToJet(3, iph, ichk);
				fTR->PhoGood[iph] += 1000;
			}
		}
		if (fTR->JGood[ichk]%10 == 0) fTR->JGood[ichk] += CleanJet(ichk, 2);
		if (fTR->JGood[ichk]%10 > 0 && fTR->JPt[ichk] > fClean_JetBadHardPtmin) fTR->GoodEvent = 4;
	}

	// Check the event
	// (Should do CleanEvent only on remaining objects)
	fTR->GoodEvent = CleanEvent();
	if( fNMuClean + fNElClean + fNJClean <= 0 || fTR->GoodEvent%10 != 0 ){
		fNMuClean = 0;
		fNElClean = 0;
		fNJClean = 0;
		return;
	}

// Check the MET
// (Should do CleanMET only on remaining objects)
// Could we choose one of the MET types?
// Or do we let the user choose???
	fTR->GoodEvent += 10   * CleanMET(fTR->TCMET, fTR->TCMETphi);
	fTR->GoodEvent += 100  * CleanMET(fTR->MuJESCorrMET, fTR->MuJESCorrMETphi);
	fTR->GoodEvent += 1000 * CleanMET(fTR->PFMET, fTR->PFMETphi);
	if( fTR->GoodEvent/10 == 111 ){
		fNMuClean = 0;
		fNElClean = 0;
		fNJClean = 0;
	}
	return;
}

void TreeCleaner::AddToJet(int ipart, int ichk, int iJet){
// adds an object (e or mu) to its nearest jet
// ipart = 1 for muon
//       = 2 for electron
//       = 3 for photon
//       = 4 for jet
// ichk = index of the object to be added
// iJet = index of the jet

	if( ichk >= 0 && iJet >= 0 ){
		double pxadd, pyadd, pzadd, eadd, ptadd, emadd;
		if( ipart == 1 ){               // muon
			pxadd = fTR->MuPx[ichk];
			pyadd = fTR->MuPy[ichk];
			pzadd = fTR->MuPz[ichk];
			eadd  = fTR->MuE[ichk];
			ptadd = fTR->MuPt[ichk];
			emadd = 0.;
		} else if( ipart == 2 ){        // electron
			pxadd = fTR->ElPx[ichk] - fTR->ElSharedPx[ichk];
			pyadd = fTR->ElPy[ichk] - fTR->ElSharedPy[ichk];
			pzadd = fTR->ElPz[ichk] - fTR->ElSharedPz[ichk];
			eadd  = fTR->ElE[ichk]  - fTR->ElSharedEnergy[ichk];
			ptadd = fTR->ElPt[ichk];
			double sharedpt = sqrt(fTR->ElSharedPx[ichk]*fTR->ElSharedPx[ichk]
				+ fTR->ElSharedPy[ichk]*fTR->ElSharedPy[ichk]);
			emadd = fTR->ElPt[ichk] - sharedpt;
			if (emadd < 0.) emadd = 0.;
		} else if( ipart == 3 ){        // photon
			pxadd = fTR->PhoPx[ichk] - fTR->PhoSharedPx[ichk];
			pyadd = fTR->PhoPy[ichk] - fTR->PhoSharedPy[ichk];
			pzadd = fTR->PhoPz[ichk] - fTR->PhoSharedPz[ichk];
			eadd  = fTR->PhoEnergy[ichk]  - fTR->PhoSharedEnergy[ichk];
			ptadd = 0.;
			double sharedpt = sqrt(fTR->PhoSharedPx[ichk]*fTR->PhoSharedPx[ichk]
				+ fTR->PhoSharedPy[ichk]*fTR->PhoSharedPy[ichk]);
			emadd = fTR->PhoPt[ichk] - sharedpt;
			if (emadd < 0.) emadd = 0.;
		} else if( ipart == 4 ){        // jet
			pxadd = fTR->JPx[ichk];
			pyadd = fTR->JPy[ichk];
			pzadd = fTR->JPz[ichk];
			eadd  = fTR->JE[ichk];
			ptadd = fTR->JChfrac[ichk] * fTR->JPt[ichk];
			emadd = fTR->JEMfrac[ichk] * fTR->JPt[ichk];
			if (emadd < 0.) emadd = 0.;
		}

		if (eadd <= 0.) return;
		fTR->JPx[iJet] += pxadd;
		fTR->JPy[iJet] += pyadd;
		fTR->JPz[iJet] += pzadd;
		fTR->JE[iJet] += eadd;
	// ??? or do we want to keep the jets massless?
		fTR->JPt[iJet] = sqrt(fTR->JPx[iJet]*fTR->JPx[iJet] + fTR->JPy[iJet]*fTR->JPy[iJet]);
		fTR->JEt[iJet] = fTR->JPt[iJet];
		if( fabs(fTR->JPz[iJet]) < 1.0e-5 ) fTR->JEta[iJet] = 0.;
		else {
			double theta = atan(fTR->JPt[iJet]/fTR->JPz[iJet]);
			if( theta < 0. ) theta += 3.141592654;
			fTR->JEta[iJet] =  -log(tan(0.5*theta));
		}
		fTR->JPhi[iJet] =  atan2(fTR->JPy[iJet], fTR->JPx[iJet]);
		if (fTR->JChfrac[iJet] >= 0.)fTR->JChfrac[iJet] += ptadd / fTR->JPt[iJet];
		fTR->JEMfrac[iJet] += emadd / fTR->JPt[iJet];
		if (fTR->JEMfrac[iJet] > 1.) fTR->JEMfrac[iJet] = 1.;
	}
	return;
}

void TreeCleaner::SubtrFromJet(int ipart, int ichk, int iJet){
// subtracts an object from its nearest jet
// ipart = 1 for muon
//       = 2 for electron
//       = 3 for photon
//       = 4 for jet
// ichk = index of the object to be subtracted
// iJet = index of the jet

  double jEcorr = fTR->JEcorr[ichk];
  if( ichk >= 0 && iJet >= 0 ){             // muon
		double pxadd, pyadd, pzadd, eadd, ptadd, emadd;
		if( ipart == 1 ){
			pxadd = fTR->MuPx[ichk];
			pyadd = fTR->MuPy[ichk];
			pzadd = fTR->MuPz[ichk];
			eadd  = fTR->MuE[ichk];
			ptadd = fTR->MuPt[ichk];
			emadd = 0.;
		} else if( ipart == 2 ){   // electron
			pxadd = fTR->ElSharedPx[ichk] * jEcorr;
			pyadd = fTR->ElSharedPy[ichk] * jEcorr;
			pzadd = fTR->ElSharedPz[ichk] * jEcorr;
			eadd  = fTR->ElSharedEnergy[ichk] * jEcorr;
			ptadd = fTR->ElPt[ichk];
			double emshpt = sqrt(fTR->ElSharedPx[ichk]*fTR->ElSharedPx[ichk]
				+ fTR->ElSharedPy[ichk]*fTR->ElSharedPy[ichk]);
			emadd = fTR->ElPt[ichk] - emshpt;
			if (emadd < 0.) emadd = 0.;
		} else if( ipart == 3 ){   // photon
			pxadd = fTR->PhoSharedPx[ichk] * jEcorr;
			pyadd = fTR->PhoSharedPy[ichk] * jEcorr;
			pzadd = fTR->PhoSharedPz[ichk] * jEcorr;
			eadd  = fTR->PhoSharedEnergy[ichk] * jEcorr;
			ptadd = 0.;
			double emshpt = sqrt(fTR->PhoSharedPx[ichk]*fTR->PhoSharedPx[ichk]
				+ fTR->PhoSharedPy[ichk]*fTR->PhoSharedPy[ichk]);
			emadd = fTR->PhoPt[ichk] - emshpt;
			if (emadd < 0.) emadd = 0.;
		} else if( ipart == 4 ){   // jet
			pxadd = fTR->JPx[ichk];
			pyadd = fTR->JPy[ichk];
			pzadd = fTR->JPz[ichk];
			eadd  = fTR->JE[ichk];
			ptadd = fTR->JChfrac[ichk] * fTR->JPt[ichk];
			emadd = fTR->JEMfrac[ichk] * fTR->JPt[ichk];
			if (emadd < 0.) emadd = 0.;
		}

		fTR->JPx[iJet] -= pxadd;
		fTR->JPy[iJet] -= pyadd;
		fTR->JPz[iJet] -= pzadd;
		fTR->JE[iJet] -= eadd;
	// ??? or do we want to keep the jets massless?
		double jmom = sqrt (fTR->JPx[iJet]*fTR->JPx[iJet] + fTR->JPy[iJet]*fTR->JPy[iJet] + fTR->JPz[iJet]*fTR->JPz[iJet]);
		if( fTR->JE[iJet] < jmom ){
			double egy = fTR->JE[iJet];
			if (egy <= 0.){ egy = 0.001;}
			double scale = egy / jmom;
			fTR->JPx[iJet] *= scale;
			fTR->JPy[iJet] *= scale;
			fTR->JPz[iJet] *= scale;
			fTR->JE[iJet] = egy;
		}
		double jPtOrigin = fTR->JPt[iJet];
		fTR->JPt[iJet] = sqrt(fTR->JPx[iJet]*fTR->JPx[iJet] + fTR->JPy[iJet]*fTR->JPy[iJet]);
		fTR->JEt[iJet] = fTR->JPt[iJet];
		if( fabs(fTR->JPz[iJet]) <1.0e-5 ) fTR->JEta[iJet] = 0.;
		else {
			double theta = atan(fTR->JPt[iJet]/fTR->JPz[iJet]);
			if( theta < 0. ) theta += 3.141592654;
			fTR->JEta[iJet] = -log(tan(0.5*theta));
		}
		fTR->JPhi[iJet] = atan2(fTR->JPy[iJet], fTR->JPx[iJet]);
		fTR->JChfrac[iJet] = (fTR->JChfrac[iJet]*jPtOrigin - ptadd) / fTR->JPt[iJet];
		fTR->JEMfrac[iJet] = (fTR->JEMfrac[iJet]*jPtOrigin - emadd) / fTR->JPt[iJet];
	}
	return;
}

int TreeCleaner::CleanEvent(void){
// To veto events from beam halo, cosmics or noise
// will also need the primary vertex
// tests on Fem and Ftrk (careful for mu-mu and e-e)
// needs to be modified when photons become available
// returns iBad = 0 event is good
//              = 1 event is empty
//              = 2 for too small EM fraction
//              = 3 for bad Trk pT fraction
// test that data still exist
	if( fNMuClean + fNElClean + fNJClean <= 0 ) return 1;

// test total Fem and Ftrk in event
	int nPhot = 0;
	int nChObj = 0;
	double pt_mu = 0.;
	double pt_track = 0.;
	double et_em = 0.;
	double et_had = 0.;
	double ptsum = 0.;
	for( int i = 0; i < fTR->NMus; ++i ){
		if(fTR->MuGood[i] != 0) continue;
		pt_mu += fTR->MuPt[i];
		pt_track += fTR->MuPt[i];
		nChObj++;
	}
	for( int i = 0; i < fTR->NEles; ++i ){
		if(fTR->ElGood[i] != 0) continue;
		pt_track += fTR->ElPt[i];
		double em = fTR->ElEt[i];
		if (em < 0.) em = 0.;
		et_em += em;
		ptsum += em;
		nChObj++;
	}
	for( int i = 0; i < fTR->NPhotons; ++i ){
		if(fTR->PhoGood[i] != 0) continue;
		double em = fTR->PhoPt[i];
		if (em < 0.) em = 0.;
		et_em += em;
		ptsum += em;
	}
	for( int i = 0; i < fTR->NJets; ++i ){
		if(fTR->JGood[i] != 0) continue;
		double pt = fTR->JChfrac[i] * fTR->JPt[i];
		if (pt < 0.) pt = 0.;
		double em = fTR->JEMfrac[i] * fTR->JEt[i];
		if (em < 0.) em = 0.;
		double had = fTR->JEt[i] - em;
		if (had < 0.) had = 0.;
		pt_track += pt;
		et_em    += em;
		et_had   += had;
		if( fTR->JChfrac[i] > 0. ) ptsum += pt;
		if( fTR->JChfrac[i] > 0. ) nChObj++;
	}

	double fracCh = 0.;
	double fracEm = 0.;
	if( et_em + et_had <= 0. ){
		if( fNMuClean < 1 ) return 1; // bad prim vtx will trigger this...
		fracCh = 1.;
		fracEm = 1.;
	} else {
		fracCh = pt_track / (ptsum + pt_mu);
		fracEm = et_em / (et_em + et_had + pt_mu);
	}
	if( fracEm < fClean_FracEmmin ) return 2;
	if( fracCh < fClean_FracChmin && (nPhot < 1 || nChObj > 0) ) return 3;
	return fTR->GoodEvent;
}

int TreeCleaner::CleanMET(double met, double metphi){
// The MET should not be aligned with any jet
// and should not be along one and opposite to the other of the 2 leading jets
// to veto QCD events with jet "mismeasurements"
// (see Jet + MET analysis in ptdr2)
// returns iBad = 0 for good MET
//              = 1 for MET aligned with a jet
//              = 2 for MET within Rij limit
	int iBad = 0;
	if( met == 0 ) return 0;
	// Care only if MET is large enough
	if (met < fClean_METmin) return 0;

	double etmax1 = 0.;
	double etmax2 = 0.;
	int imax1 = -1;
	int imax2 = -1;

	// Loop over all jets
	for( int i = 0; i < fTR->NJets; ++i ){
		// Reject if the MET is along the jet
		if(fTR->JGood[i] != 0) continue;
		double dPhi =  Util::DeltaPhi(fTR->JPhi[i], metphi);
		if( dPhi < fClean_dPhiJetMETmin ) return 1;
		// Else, pick up the 2 leading jets
		if( fTR->JPt[i] > etmax1 ){
			etmax2 = etmax1;
			imax2  = imax1;
			etmax1 = fTR->JPt[i];
			imax1  = i;
		} else if( fTR->JPt[i] > etmax2 ){
			etmax2 = fTR->JPt[i];
			imax2  = i;
		}
	}

	// Check dR12 and dR21
	if( imax2 >= 0 ){
		double dPhi1 = Util::DeltaPhi(fTR->JPhi[imax1], metphi );
		double dPhi2 = Util::DeltaPhi(fTR->JPhi[imax2], metphi );
		double pi = 3.141592654;
		fR12 = sqrt(dPhi1*dPhi1 + (pi-dPhi2)*(pi-dPhi2) );
		fR21 = sqrt(dPhi2*dPhi2 + (pi-dPhi1)*(pi-dPhi1) );
		if( fR12 < fClean_dR12min || fR21 < fClean_dR21min ) return 2;
	}
	return 0;
}

int TreeCleaner::FindNearestJet(double eta, double phi){
// Looks for the nearest jet in deltaR to a given object (e or mu or photon)
// and returns its index
// returns -1 if no nearest jet

	int iJetMin = -1;

	double deltaRmin = 999.;
	for(int i = 0; i < fTR->NJets; ++i){
		double deltaR = Util::GetDeltaR(eta, fTR->JEta[i], phi, fTR->JPhi[i]);
		if (deltaR < deltaRmin){
			deltaRmin = deltaR;
			iJetMin = i;
		}
	}
	return iJetMin;
}

void TreeCleaner::StatInit(const char* filename){
// initializes the cleaning statistics
// to be called once at the beginning of the job
	fNumTotEvt               = 0;
	fNumTotEvtReject         = 0;
	fNumTotEvtEmpty          = 0;
	fNumTotEvtCleanEmpty     = 0;
	fNumTotEvtLtFem          = 0;
	fNumTotEvtLtFch          = 0;
	fNumTotEvtPfMETJet       = 0;
	fNumTotEvtPfMETRij       = 0;
	fNumTotEvtCaMETJet       = 0;
	fNumTotEvtCaMETRij       = 0;
	fNumTotEvtTcMETJet       = 0;
	fNumTotEvtTcMETRij       = 0;
	fNumTotEvtBadHardJet     = 0;

	fNumTotMuons             = 0;  
	fNumTotMuonGoodIso       = 0;  
	fNumTotMuonGoodNonIso    = 0;  
	fNumTotMuonBadIso        = 0;  
	fNumTotMuonBadNonIso     = 0;  
	fNumTotMuonDupl          = 0;
	fNumTotMuonNotPrimaryTrk = 0;
	fNumTotMuonNotClean      = 0;
	fNumTotMuonBadDpop       = 0;
	fNumTotMuonBadChi2       = 0;  
	fNumTotMuonBadNhit       = 0;  
	fNumTotMuonBadGlTr       = 0;  

	fNumTotElectrons         = 0;
	fNumTotElecGoodIso       = 0;
	fNumTotElecGoodNonIso    = 0;
	fNumTotElecBadIso        = 0;
	fNumTotElecBadNonIso     = 0;
	fNumTotElecDupl          = 0;
	fNumTotElecNotPrimaryTrk = 0;
	fNumTotElecNotClean      = 0;
	fNumTotElecBadHoE        = 0;
	fNumTotElecBadShsh       = 0;
	fNumTotElecBadTmat       = 0;
	fNumTotElecBadSpik       = 0;
	fNumTotElecBadHits       = 0;
	fNumTotElecBadConv       = 0;

	fNumTotPhotons           = 0;
	fNumTotPhotGoodIso       = 0;
	fNumTotPhotGoodNonIso    = 0;
	fNumTotPhotBadIso        = 0;
	fNumTotPhotBadNonIso     = 0;
	fNumTotPhotDupl          = 0;
	fNumTotPhotNotClean      = 0;
	fNumTotPhotBadHoE        = 0;
	fNumTotPhotBadShsh       = 0;
	fNumTotPhotBadSpik       = 0;

	fNumTotJets              = 0;  
	fNumTotJetGood           = 0;  
	fNumTotJetBad            = 0;  
	fNumTotJetDuplElJet      = 0;
	fNumTotJetDuplPhoJet     = 0;
	fNumTotJetNotPrimaryTrk  = 0;
	fNumTotJetNotClean       = 0;
	fNumTotJetPgtE           = 0;
	fNumTotJetGtFem          = 0;
	fNumTotJetLtFem          = 0;
	fNumTotJetLtFch          = 0;
	fNumTotJetLtn90hits      = 0;
	fNumTotJetGtfHPD         = 0;
	fNumTotJetGtfRBX         = 0;
	fNumTotBJets             = 0;  

	// initialize cleaning statistics histogram
	fHstatFile = new TFile(fOutputDir + TString(filename), "RECREATE");
	fHstatHistos  = new TH1D("CleanStats", "cleaning statistics", 100, 0, 100);

	return;
}

void TreeCleaner::StatFill(void){
// accumulates the cleaning statistics
// to be called once per event, after cleaning is done

	fNumTotEvt++;
	int nObjects = 0;

	fNumTotMuons += fTR->NMus;
	for( int i=0; i < fTR->NMus; ++i ){
		if(      fTR->MuGood[i] == 0 && fTR->MuIsIso[i] == 1 ){ fNumTotMuonGoodIso++; nObjects++; }
		else if( fTR->MuGood[i] == 0 && fTR->MuIsIso[i] == 0 ) fNumTotMuonGoodNonIso++;
		else if( fTR->MuGood[i] != 0 && fTR->MuIsIso[i] == 1 ) fNumTotMuonBadIso++;
		else if( fTR->MuGood[i] != 0 && fTR->MuIsIso[i] == 0 ) fNumTotMuonBadNonIso++;
		if( fTR->MuGood[i] !=0 ){
			if( fTR->MuGood[i] / 100 != 0 ) fNumTotMuonDupl++;
			if( fTR->MuGood[i] % 100 / 10 != 0 ) fNumTotMuonNotPrimaryTrk++;
			int muClean = fTR->MuGood[i] % 10;
			if( muClean != 0 ){
				fNumTotMuonNotClean++;
				if( muClean == 1 ) fNumTotMuonBadDpop++;
				if( muClean == 2 ) fNumTotMuonBadChi2++;
				if( muClean == 3 ) fNumTotMuonBadNhit++;
				if( muClean == 4 ) fNumTotMuonBadGlTr++;
			}
		}
	}

	fNumTotElectrons += fTR->NEles;
	for( int i=0; i < fTR->NEles; ++i ){
		if(      fTR->ElGood[i] == 0 && fTR->ElIsIso[i] == 1 ){ fNumTotElecGoodIso++; nObjects++; }
		else if( fTR->ElGood[i] == 0 && fTR->ElIsIso[i] == 0 ) fNumTotElecGoodNonIso++;
		else if( fTR->ElGood[i] != 0 && fTR->ElIsIso[i] == 1 ) fNumTotElecBadIso++;
		else if( fTR->ElGood[i] != 0 && fTR->ElIsIso[i] == 0 ) fNumTotElecBadNonIso++;
		if( fTR->ElGood[i] != 0 ){
			if( fTR->ElGood[i] / 100 != 0 ) fNumTotElecDupl++;
			if( fTR->ElGood[i] % 100 / 10 != 0 ) fNumTotElecNotPrimaryTrk++;
			int elClean = fTR->ElGood[i] % 10;
			if( elClean != 0 ){
				fNumTotElecNotClean++;
				if( elClean == 1 ) fNumTotElecBadHoE++;
				if( elClean == 2 ) fNumTotElecBadShsh++;
				if( elClean == 3 ) fNumTotElecBadTmat++;
				if( elClean == 4 ) fNumTotElecBadSpik++;
				if( elClean == 5 ) fNumTotElecBadHits++;
				if( elClean == 6 ) fNumTotElecBadConv++;
			}
		}
	}

	fNumTotPhotons += fTR->NPhotons;
	for( int i=0; i < fTR->NPhotons; ++i ){
		if(      fTR->PhoGood[i] == 0 && fTR->PhoIsIso[i] == 1 ){ fNumTotPhotGoodIso++; nObjects++; }
		else if( fTR->PhoGood[i] == 0 && fTR->PhoIsIso[i] == 0 ) fNumTotPhotGoodNonIso++;
		else if( fTR->PhoGood[i] != 0 && fTR->PhoIsIso[i] == 1 ) fNumTotPhotBadIso++;
		else if( fTR->PhoGood[i] != 0 && fTR->PhoIsIso[i] == 0 ) fNumTotPhotBadNonIso++;
		if( fTR->PhoGood[i] != 0 ){
			int phClean = fTR->PhoGood[i] % 10;
			if( phClean != 0 ){
				fNumTotPhotNotClean++;
				if( phClean == 1 ) fNumTotPhotBadHoE++;
				if( phClean == 2 ) fNumTotPhotBadShsh++;
				if( phClean == 3 ) fNumTotPhotBadSpik++;
			}
		}
	}

	fNumTotJets += fTR->NJets;
	for( int i=0; i < fTR->NJets; ++i ){
		if( fTR->JGood[i] == 0 ){ fNumTotJetGood++; nObjects++; }
		else {
			fNumTotJetBad++;
			if( fTR->JGood[i] / 100 == 1 ) fNumTotJetDuplElJet++;
			if( fTR->JGood[i] / 100 == 2 ) fNumTotJetDuplPhoJet++;
			if( fTR->JGood[i] % 100 / 10 != 0 ) fNumTotJetNotPrimaryTrk++;
			int jetClean = fTR->JGood[i] % 10;
			if( jetClean != 0 ){
				fNumTotJetNotClean++;
				if( jetClean == 1 ) fNumTotJetPgtE++;
				if( jetClean == 2 ) fNumTotJetGtFem++;
				if( jetClean == 3 ) fNumTotJetLtFem++;
				if( jetClean == 4 ) fNumTotJetLtFch++;
				if( jetClean == 5 ) fNumTotJetLtn90hits++;
				if( jetClean == 6 ) fNumTotJetGtfHPD++;
				if( jetClean == 7 ) fNumTotJetGtfRBX++;
			}
		}
	}

	if( nObjects <= 0) fNumTotEvtEmpty++;
	int evClean = fTR->GoodEvent % 10;
	if( evClean != 0 || nObjects <= 0 ) fNumTotEvtReject++;
	if( evClean == 1 ) fNumTotEvtCleanEmpty++;
	if( evClean == 2 ) fNumTotEvtLtFem++;
	if( evClean == 3 ) fNumTotEvtLtFch++;
	if( evClean == 4 ) fNumTotEvtBadHardJet++;

	evClean = fTR->GoodEvent / 1000;
	if( evClean == 1 ) fNumTotEvtPfMETJet++;
	if( evClean == 2 ) fNumTotEvtPfMETRij++;
	evClean = fTR->GoodEvent % 1000 / 100;
	if( evClean == 1 ) fNumTotEvtCaMETJet++;
	if( evClean == 2 ) fNumTotEvtCaMETRij++;
	evClean = fTR->GoodEvent % 100 / 10;
	if( evClean == 1 ) fNumTotEvtTcMETJet++;
	if( evClean == 2 ) fNumTotEvtTcMETRij++;
	return;
}

void TreeCleaner::StatPrint(void){
// prints the cleaning statistics
// to be called once at the end of the job


	cout << endl;
	cout << "Cleaning statistics from TreeCleaner " << endl;
	// cout << "  " << fNumTotElecGoodIso << "  " << fNumTotElecGoodNonIso
	//         << "  " << fNumTotElecBadIso  << "  " << fNumTotElecBadNonIso << endl;

	cout << endl;
	cout << " Total number of events processed = " << fNumTotEvt << endl;

// Statistics for events
	if (fNumTotEvt > 0) {
		cout << endl;
		cout << "   events accepted                  = " 
			<< fNumTotEvt-fNumTotEvtReject << endl;
		cout << "   events rejected (total)          = " << fNumTotEvtReject
			<< "  = " << 100.*(float)fNumTotEvtReject / (float)fNumTotEvt << " %" << endl;
		cout << endl;
		cout << "    empty after cleaning+Iso        = " << fNumTotEvtEmpty
			<< "  = " << 100.*(float)fNumTotEvtEmpty / (float)fNumTotEvt << " %" << endl;
		cout << "    empty after cleaning            = " << fNumTotEvtCleanEmpty
			<< "  = " << 100.*(float)fNumTotEvtCleanEmpty / (float)fNumTotEvt << " %" << endl;
		cout << "    too small em fraction           = " << fNumTotEvtLtFem
			<< "  = " << 100.*(float)fNumTotEvtLtFem / (float)fNumTotEvt << " %" << endl;
		cout << "    too small track pT fraction     = " << fNumTotEvtLtFch
			<< "  = " << 100.*(float)fNumTotEvtLtFch / (float)fNumTotEvt << " %" << endl;
		cout << "    jet with pT > " << fClean_JetBadHardPtmin <<  " and bad        = " << fNumTotEvtBadHardJet
			<< "  = " << 100.*(float)fNumTotEvtBadHardJet / (float)fNumTotEvt << " %" << endl;
		cout << endl;
		cout << "    caloMET aligned with jet        = " << fNumTotEvtCaMETJet
			<< "  = " << 100.*(float)fNumTotEvtCaMETJet  / (float)fNumTotEvt << " %" << endl;
		cout << "    caloMET bad Rij                 = " << fNumTotEvtCaMETRij
			<< "  = " << 100.*(float)fNumTotEvtCaMETRij  / (float)fNumTotEvt << " %" << endl;
		cout << "    tcMET aligned with jet          = " << fNumTotEvtCaMETJet
			<< "  = " << 100.*(float)fNumTotEvtCaMETJet  / (float)fNumTotEvt << " %" << endl;
		cout << "    tcMET bad Rij                   = " << fNumTotEvtCaMETRij
			<< "  = " << 100.*(float)fNumTotEvtCaMETRij  / (float)fNumTotEvt << " %" << endl;
		cout << "    pfMET aligned with jet          = " << fNumTotEvtCaMETJet
			<< "  = " << 100.*(float)fNumTotEvtCaMETJet  / (float)fNumTotEvt << " %" << endl;
		cout << "    pfMET bad Rij                   = " << fNumTotEvtCaMETRij
			<< "  = " << 100.*(float)fNumTotEvtCaMETRij  / (float)fNumTotEvt << " %" << endl ;
	}

// Statistics for muons
	cout << endl;
	cout << " Total number of muons              = " << fNumTotMuons << endl;
	if (fNumTotMuons > 0) {
		cout << "  muons good+Iso                    = " << fNumTotMuonGoodIso
			<< "  = " << 100.*(float)fNumTotMuonGoodIso / (float)fNumTotMuons << " %" << endl;
		cout << "  muons good+non-Iso                = " << fNumTotMuonGoodNonIso
			<< "  = " << 100.*(float)fNumTotMuonGoodNonIso / (float)fNumTotMuons << " %" << endl;
		cout << "  muons bad+Iso                     = " << fNumTotMuonBadIso
			<< "  = " << 100.*(float)fNumTotMuonBadIso / (float)fNumTotMuons << " %" << endl;
		cout << "  muons bad+non-Iso                 = " << fNumTotMuonBadNonIso
			<< "  = " << 100.*(float)fNumTotMuonBadNonIso / (float)fNumTotMuons << " %" << endl;
		cout << endl;
		int mubad = fNumTotMuonBadIso + fNumTotMuonBadNonIso;
		cout << "  muons bad total                   = " << mubad
			<< "  = " << 100.*(float)mubad / (float)fNumTotMuons << " %" << endl;
		cout << "   muons duplicated                 = " << fNumTotMuonDupl
			<< "  = " << 100.*(float)fNumTotMuonDupl / (float)fNumTotMuons << " %" << endl;
		cout << "   muons not from Primary Vertex    = " << fNumTotMuonNotPrimaryTrk
			<< "  = " << 100.*(float)fNumTotMuonNotPrimaryTrk / (float)fNumTotMuons << " %" << endl;
		cout << "   muons not clean                  = " << fNumTotMuonNotClean
			<< "  = " << 100.*(float)fNumTotMuonNotClean / (float)fNumTotMuons << " %" << endl;
		cout << "    muons with bad dP/P             = " << fNumTotMuonBadDpop
			<< "  = " << 100.*(float)fNumTotMuonBadDpop / (float)fNumTotMuons << " %" << endl;
		cout << "    muons with bad chisquared       = " << fNumTotMuonBadChi2
			<< "  = " << 100.*(float)fNumTotMuonBadChi2 / (float)fNumTotMuons << " %" << endl;
		cout << "    muons with bad #tracker hits    = " << fNumTotMuonBadNhit
			<< "  = " << 100.*(float)fNumTotMuonBadNhit / (float)fNumTotMuons << " %" << endl;
		cout << "    muons not global+tracker        = " << fNumTotMuonBadGlTr
			<< "  = " << 100.*(float)fNumTotMuonBadGlTr / (float)fNumTotMuons << " %" << endl;
	}

// Statistics for electrons
	cout << endl;
	cout << " Total number of electrons          = " << fNumTotElectrons << endl;
	if (fNumTotElectrons > 0) {
		cout << "  elecs good+Iso                    = " << fNumTotElecGoodIso
			<< "  = " << 100.*(float)fNumTotElecGoodIso / (float)fNumTotElectrons << " %" << endl;
		cout << "  elecs good+non-Iso                = " << fNumTotElecGoodNonIso
			<< "  = " << 100.*(float)fNumTotElecGoodNonIso / (float)fNumTotElectrons << " %" << endl;
		cout << "  elecs bad+Iso                     = " << fNumTotElecBadIso
			<< "  = " << 100.*(float)fNumTotElecBadIso / (float)fNumTotElectrons << " %" << endl;
		cout << "  elecs bad+non-Iso                 = " << fNumTotElecBadNonIso
			<< "  = " << 100.*(float)fNumTotElecBadNonIso / (float)fNumTotElectrons << " %" << endl;
		cout << endl;
		int elebad = fNumTotElecBadIso + fNumTotElecBadNonIso;
		cout << "  elecs bad total                   = " << elebad
			<< "  = " << 100.*(float)elebad / (float)fNumTotElectrons << " %" << endl;
		cout << "   elecs duplicated                 = " << fNumTotElecDupl
			<< "  = " << 100.*(float)fNumTotElecDupl / (float)fNumTotElectrons << " %" << endl;
		cout << "   elecs not from Primary Vertex    = " << fNumTotElecNotPrimaryTrk
			<< "  = " << 100.*(float)fNumTotElecNotPrimaryTrk / (float)fNumTotElectrons << " %" << endl;
		cout << "   elecs not clean                  = " << fNumTotElecNotClean
			<< "  = " << 100.*(float)fNumTotElecNotClean / (float)fNumTotElectrons << " %" << endl;
		cout << "    elecs with bad H/E              = " << fNumTotElecBadHoE
			<< "  = " << 100.*(float)fNumTotElecBadHoE / (float)fNumTotElectrons << " %" << endl;
		cout << "    elecs with bad shower shape     = " << fNumTotElecBadShsh
			<< "  = " << 100.*(float)fNumTotElecBadShsh / (float)fNumTotElectrons << " %" << endl;
		cout << "    elecs with bad track matching   = " << fNumTotElecBadTmat
			<< "  = " << 100.*(float)fNumTotElecBadTmat / (float)fNumTotElectrons << " %" << endl;
		cout << "    elecs from spike                = " << fNumTotElecBadSpik
			<< "  = " << 100.*(float)fNumTotElecBadSpik / (float)fNumTotElectrons << " %" << endl;
		cout << "    elecs with missing inner hits   = " << fNumTotElecBadHits
			<< "  = " << 100.*(float)fNumTotElecBadHits / (float)fNumTotElectrons << " %" << endl;
		cout << "    elecs with conversion track     = " << fNumTotElecBadConv
			<< "  = " << 100.*(float)fNumTotElecBadConv / (float)fNumTotElectrons << " %" << endl;
	}

// Statistics for photons
	cout << endl;
	cout << " Total number of photons            = " << fNumTotPhotons << endl;
	if (fNumTotPhotons > 0) {
		cout << "  photons good+Iso                  = " << fNumTotPhotGoodIso
			<< "  = " << 100.*(float)fNumTotPhotGoodIso / (float)fNumTotPhotons << " %" << endl;
		cout << "  photons good+non-Iso              = " << fNumTotPhotGoodNonIso
			<< "  = " << 100.*(float)fNumTotPhotGoodNonIso / (float)fNumTotPhotons << " %" << endl;
		cout << "  photons bad+Iso                   = " << fNumTotPhotBadIso
			<< "  = " << 100.*(float)fNumTotPhotBadIso / (float)fNumTotPhotons << " %" << endl;
		cout << "  photons bad+non-Iso               = " << fNumTotPhotBadNonIso
			<< "  = " << 100.*(float)fNumTotPhotBadNonIso / (float)fNumTotPhotons << " %" << endl;
		cout << endl;
		int photbad = fNumTotPhotBadIso + fNumTotPhotBadNonIso;
		cout << "  photons bad total                 = " << photbad
			<< "  = " << 100.*(float)photbad / (float)fNumTotPhotons << " %" << endl;
		cout << "   photons duplicated               = " << fNumTotPhotDupl
			<< "  = " << 100.*(float)fNumTotPhotDupl / (float)fNumTotPhotons << " %" << endl;
		cout << "   photons not clean                = " << fNumTotPhotNotClean
			<< "  = " << 100.*(float)fNumTotPhotNotClean / (float)fNumTotPhotons << " %" << endl;
		cout << "    photons with bad H/E            = " << fNumTotPhotBadHoE
			<< "  = " << 100.*(float)fNumTotPhotBadHoE / (float)fNumTotPhotons << " %" << endl;
		cout << "    photons with bad shower shape   = " << fNumTotPhotBadShsh
			<< "  = " << 100.*(float)fNumTotPhotBadShsh / (float)fNumTotPhotons << " %" << endl;
		cout << "    photons from spike              = " << fNumTotPhotBadSpik
			<< "  = " << 100.*(float)fNumTotPhotBadSpik / (float)fNumTotPhotons << " %" << endl;
	}

// Statistics for jets
	cout << endl;
	cout << " Total number of jets               = " << fNumTotJets << endl;
	if (fNumTotJets > 0) {
		cout << "  jets good                         = " << fNumTotJetGood
			<< "  = " << 100.*(float)fNumTotJetGood / (float)fNumTotJets << " %" << endl;
		cout << "  jets bad                          = " << fNumTotJetBad
			<< "  = " << 100.*(float)fNumTotJetBad / (float)fNumTotJets << " %" << endl;
		cout << endl;
		cout << "   jets duplicated with elec        = " << fNumTotJetDuplElJet
			<< "  = " << 100.*(float)fNumTotJetDuplElJet / (float)fNumTotJets << " %" << endl;
		cout << "   jets duplicated with photon      = " << fNumTotJetDuplPhoJet
			<< "  = " << 100.*(float)fNumTotJetDuplPhoJet / (float)fNumTotJets << " %" << endl;
		cout << "   jets not from Primary Vertex     = " << fNumTotJetNotPrimaryTrk
			<< "  = " << 100.*(float)fNumTotJetNotPrimaryTrk / (float)fNumTotJets << " %" << endl;
		cout << "   jets not clean                   = " << fNumTotJetNotClean
			<< "  = " << 100.*(float)fNumTotJetNotClean / (float)fNumTotJets << " %" << endl;
		cout << "    jets with P > E                 = " << fNumTotJetPgtE
			<< "  = " << 100.*(float)fNumTotJetPgtE / (float)fNumTotJets << " %" << endl;
		cout << "    jets with too large Fem         = " << fNumTotJetGtFem
			<< "  = " << 100.*(float)fNumTotJetGtFem / (float)fNumTotJets << " %" << endl;
		cout << "    jets with too small Fem         = " << fNumTotJetLtFem
			<< "  = " << 100.*(float)fNumTotJetLtFem / (float)fNumTotJets << " %" << endl;
		cout << "    jets with too small Fch         = " << fNumTotJetLtFch
			<< "  = " << 100.*(float)fNumTotJetLtFch / (float)fNumTotJets << " %" << endl;
		cout << "    jets with too small n90Hits     = " << fNumTotJetLtn90hits
			<< "  = " << 100.*(float)fNumTotJetLtn90hits / (float)fNumTotJets << " %" << endl;
		cout << "    jets with too large fHPD        = " << fNumTotJetGtfHPD
			<< "  = " << 100.*(float)fNumTotJetGtfHPD / (float)fNumTotJets << " %" << endl;
		cout << "    jets with too large fRBX        = " << fNumTotJetGtfRBX
			<< "  = " << 100.*(float)fNumTotJetGtfRBX / (float)fNumTotJets << " %" << endl;
		cout << "  Total number of b-jets            = " << fNumTotBJets
			<< "  = " << 100.*(float)fNumTotBJets / (float)fNumTotJets << " %" << endl;
	}
	return;  
}

void TreeCleaner::StatWrite(TString filename){
// writes the cleaning statistics to file
// to be called once at the end of the job

	ofstream file;
	file.open(filename, ios::out);

	file << endl;
	file << "Cleaning statistics from TreeCleaner " << endl;
	// file << "  " << fNumTotElecGoodIso << "  " << fNumTotElecGoodNonIso
	//         << "  " << fNumTotElecBadIso  << "  " << fNumTotElecBadNonIso << endl;

	file << endl;
	file << " Total number of events processed = " << fNumTotEvt << endl;

// Statistics for events
	if (fNumTotEvt > 0) {
		file << endl;
		file << "   events accepted                  = " 
			<< fNumTotEvt-fNumTotEvtReject << endl;
		file << "   events rejected (total)          = " << fNumTotEvtReject
			<< "  = " << 100.*(float)fNumTotEvtReject / (float)fNumTotEvt << " %" << endl;
		file << endl;
		file << "    empty after cleaning+Iso        = " << fNumTotEvtEmpty
			<< "  = " << 100.*(float)fNumTotEvtEmpty / (float)fNumTotEvt << " %" << endl;
		file << "    empty after cleaning            = " << fNumTotEvtCleanEmpty
			<< "  = " << 100.*(float)fNumTotEvtCleanEmpty / (float)fNumTotEvt << " %" << endl;
		file << "    too small em fraction           = " << fNumTotEvtLtFem
			<< "  = " << 100.*(float)fNumTotEvtLtFem / (float)fNumTotEvt << " %" << endl;
		file << "    too small track pT fraction     = " << fNumTotEvtLtFch
			<< "  = " << 100.*(float)fNumTotEvtLtFch / (float)fNumTotEvt << " %" << endl;
		file << "    jet with pT > " << fClean_JetBadHardPtmin <<  " and bad        = " << fNumTotEvtBadHardJet
			<< "  = " << 100.*(float)fNumTotEvtBadHardJet / (float)fNumTotEvt << " %" << endl;
		file << endl;
		file << "    caloMET aligned with jet        = " << fNumTotEvtCaMETJet
			<< "  = " << 100.*(float)fNumTotEvtCaMETJet  / (float)fNumTotEvt << " %" << endl;
		file << "    caloMET bad Rij                 = " << fNumTotEvtCaMETRij
			<< "  = " << 100.*(float)fNumTotEvtCaMETRij  / (float)fNumTotEvt << " %" << endl;
		file << "    tcMET aligned with jet          = " << fNumTotEvtCaMETJet
			<< "  = " << 100.*(float)fNumTotEvtCaMETJet  / (float)fNumTotEvt << " %" << endl;
		file << "    tcMET bad Rij                   = " << fNumTotEvtCaMETRij
			<< "  = " << 100.*(float)fNumTotEvtCaMETRij  / (float)fNumTotEvt << " %" << endl;
		file << "    pfMET aligned with jet          = " << fNumTotEvtCaMETJet
			<< "  = " << 100.*(float)fNumTotEvtCaMETJet  / (float)fNumTotEvt << " %" << endl;
		file << "    pfMET bad Rij                   = " << fNumTotEvtCaMETRij
			<< "  = " << 100.*(float)fNumTotEvtCaMETRij  / (float)fNumTotEvt << " %" << endl ;
	}

// Statistics for muons
	file << endl;
	file << " Total number of muons              = " << fNumTotMuons << endl;
	if (fNumTotMuons > 0) {
		file << "  muons good+Iso                    = " << fNumTotMuonGoodIso
			<< "  = " << 100.*(float)fNumTotMuonGoodIso / (float)fNumTotMuons << " %" << endl;
		file << "  muons good+non-Iso                = " << fNumTotMuonGoodNonIso
			<< "  = " << 100.*(float)fNumTotMuonGoodNonIso / (float)fNumTotMuons << " %" << endl;
		file << "  muons bad+Iso                     = " << fNumTotMuonBadIso
			<< "  = " << 100.*(float)fNumTotMuonBadIso / (float)fNumTotMuons << " %" << endl;
		file << "  muons bad+non-Iso                 = " << fNumTotMuonBadNonIso
			<< "  = " << 100.*(float)fNumTotMuonBadNonIso / (float)fNumTotMuons << " %" << endl;
		file << endl;
		int mubad = fNumTotMuonBadIso + fNumTotMuonBadNonIso;
		file << "  muons bad total                   = " << mubad
			<< "  = " << 100.*(float)mubad / (float)fNumTotMuons << " %" << endl;
		file << "   muons duplicated                 = " << fNumTotMuonDupl
			<< "  = " << 100.*(float)fNumTotMuonDupl / (float)fNumTotMuons << " %" << endl;
		file << "   muons not from Primary Vertex    = " << fNumTotMuonNotPrimaryTrk
			<< "  = " << 100.*(float)fNumTotMuonNotPrimaryTrk / (float)fNumTotMuons << " %" << endl;
		file << "   muons not clean                  = " << fNumTotMuonNotClean
			<< "  = " << 100.*(float)fNumTotMuonNotClean / (float)fNumTotMuons << " %" << endl;
		file << "    muons with bad dP/P             = " << fNumTotMuonBadDpop
			<< "  = " << 100.*(float)fNumTotMuonBadDpop / (float)fNumTotMuons << " %" << endl;
		file << "    muons with bad chisquared       = " << fNumTotMuonBadChi2
			<< "  = " << 100.*(float)fNumTotMuonBadChi2 / (float)fNumTotMuons << " %" << endl;
		file << "    muons with bad #tracker hits    = " << fNumTotMuonBadNhit
			<< "  = " << 100.*(float)fNumTotMuonBadNhit / (float)fNumTotMuons << " %" << endl;
	}

// Statistics for electrons
	file << endl;
	file << " Total number of electrons          = " << fNumTotElectrons << endl;
	if (fNumTotElectrons > 0) {
		file << "  elecs good+Iso                    = " << fNumTotElecGoodIso
			<< "  = " << 100.*(float)fNumTotElecGoodIso / (float)fNumTotElectrons << " %" << endl;
		file << "  elecs good+non-Iso                = " << fNumTotElecGoodNonIso
			<< "  = " << 100.*(float)fNumTotElecGoodNonIso / (float)fNumTotElectrons << " %" << endl;
		file << "  elecs bad+Iso                     = " << fNumTotElecBadIso
			<< "  = " << 100.*(float)fNumTotElecBadIso / (float)fNumTotElectrons << " %" << endl;
		file << "  elecs bad+non-Iso                 = " << fNumTotElecBadNonIso
			<< "  = " << 100.*(float)fNumTotElecBadNonIso / (float)fNumTotElectrons << " %" << endl;
		file << endl;
		int elebad = fNumTotElecBadIso + fNumTotElecBadNonIso;
		file << "  elecs bad total                   = " << elebad
			<< "  = " << 100.*(float)elebad / (float)fNumTotElectrons << " %" << endl;
		file << "   elecs duplicated                 = " << fNumTotElecDupl
			<< "  = " << 100.*(float)fNumTotElecDupl / (float)fNumTotElectrons << " %" << endl;
		file << "   elecs not from Primary Vertex    = " << fNumTotElecNotPrimaryTrk
			<< "  = " << 100.*(float)fNumTotElecNotPrimaryTrk / (float)fNumTotElectrons << " %" << endl;
		file << "   elecs not clean                  = " << fNumTotElecNotClean
			<< "  = " << 100.*(float)fNumTotElecNotClean / (float)fNumTotElectrons << " %" << endl;
		file << "    elecs with bad H/E              = " << fNumTotElecBadHoE
			<< "  = " << 100.*(float)fNumTotElecBadHoE / (float)fNumTotElectrons << " %" << endl;
		file << "    elecs with bad shower shape     = " << fNumTotElecBadShsh
			<< "  = " << 100.*(float)fNumTotElecBadShsh / (float)fNumTotElectrons << " %" << endl;
		file << "    elecs with bad track matching   = " << fNumTotElecBadTmat
			<< "  = " << 100.*(float)fNumTotElecBadTmat / (float)fNumTotElectrons << " %" << endl;
	}

// Statistics for photons
	file << endl;
	file << " Total number of photons            = " << fNumTotPhotons << endl;
	if (fNumTotPhotons > 0) {
		file << "  photons good+Iso                  = " << fNumTotPhotGoodIso
			<< "  = " << 100.*(float)fNumTotPhotGoodIso / (float)fNumTotPhotons << " %" << endl;
		file << "  photons good+non-Iso              = " << fNumTotPhotGoodNonIso
			<< "  = " << 100.*(float)fNumTotPhotGoodNonIso / (float)fNumTotPhotons << " %" << endl;
		file << "  photons bad+Iso                   = " << fNumTotPhotBadIso
			<< "  = " << 100.*(float)fNumTotPhotBadIso / (float)fNumTotPhotons << " %" << endl;
		file << "  photons bad+non-Iso               = " << fNumTotPhotBadNonIso
			<< "  = " << 100.*(float)fNumTotPhotBadNonIso / (float)fNumTotPhotons << " %" << endl;
		file << endl;
		int photbad = fNumTotPhotBadIso + fNumTotPhotBadNonIso;
		file << "  photons bad total                 = " << photbad
			<< "  = " << 100.*(float)photbad / (float)fNumTotPhotons << " %" << endl;
		file << "   photons duplicated               = " << fNumTotPhotDupl
			<< "  = " << 100.*(float)fNumTotPhotDupl / (float)fNumTotPhotons << " %" << endl;
		file << "   photons not clean                = " << fNumTotPhotNotClean
			<< "  = " << 100.*(float)fNumTotPhotNotClean / (float)fNumTotPhotons << " %" << endl;
		file << "    photons with bad H/E            = " << fNumTotPhotBadHoE
			<< "  = " << 100.*(float)fNumTotPhotBadHoE / (float)fNumTotPhotons << " %" << endl;
		file << "    photons with bad shower shape   = " << fNumTotPhotBadShsh
			<< "  = " << 100.*(float)fNumTotPhotBadShsh / (float)fNumTotPhotons << " %" << endl;
	}

// Statistics for jets
	file << endl;
	file << " Total number of jets               = " << fNumTotJets << endl;
	if (fNumTotJets > 0) {
		file << "  jets good                         = " << fNumTotJetGood
			<< "  = " << 100.*(float)fNumTotJetGood / (float)fNumTotJets << " %" << endl;
		file << "  jets bad                          = " << fNumTotJetBad
			<< "  = " << 100.*(float)fNumTotJetBad / (float)fNumTotJets << " %" << endl;
		file << endl;
		file << "   jets duplicated with elec        = " << fNumTotJetDuplElJet
			<< "  = " << 100.*(float)fNumTotJetDuplElJet / (float)fNumTotJets << " %" << endl;
		file << "   jets duplicated with photon      = " << fNumTotJetDuplPhoJet
			<< "  = " << 100.*(float)fNumTotJetDuplPhoJet / (float)fNumTotJets << " %" << endl;
		file << "   jets not from Primary Vertex     = " << fNumTotJetNotPrimaryTrk
			<< "  = " << 100.*(float)fNumTotJetNotPrimaryTrk / (float)fNumTotJets << " %" << endl;
		file << "   jets not clean                   = " << fNumTotJetNotClean
			<< "  = " << 100.*(float)fNumTotJetNotClean / (float)fNumTotJets << " %" << endl;
		file << "    jets with P > E                 = " << fNumTotJetPgtE
			<< "  = " << 100.*(float)fNumTotJetPgtE / (float)fNumTotJets << " %" << endl;
		file << "    jets with too large Fem         = " << fNumTotJetGtFem
			<< "  = " << 100.*(float)fNumTotJetGtFem / (float)fNumTotJets << " %" << endl;
		file << "    jets with too small Fem         = " << fNumTotJetLtFem
			<< "  = " << 100.*(float)fNumTotJetLtFem / (float)fNumTotJets << " %" << endl;
		file << "    jets with too small Fch         = " << fNumTotJetLtFch
			<< "  = " << 100.*(float)fNumTotJetLtFch / (float)fNumTotJets << " %" << endl;
		file << "  Total number of b-jets            = " << fNumTotBJets
			<< "  = " << 100.*(float)fNumTotBJets / (float)fNumTotJets << " %" << endl;
	}

	file.close();

	return;  
}

void TreeCleaner::StatHistos(void){
// saves the cleaning statistics into histograms (from Kostas)
// to be called once at the end of the job

	fHstatHistos->SetFillColor(kRed);
	fHstatHistos->SetDrawOption("hbar");
	fHstatHistos->GetXaxis()->SetBinLabel(1,  "TotEvts");           fHstatHistos->SetBinContent(1,  fNumTotEvt);
	fHstatHistos->GetXaxis()->SetBinLabel(2,  "EvtReject");         fHstatHistos->SetBinContent(2,  fNumTotEvtReject);
	fHstatHistos->GetXaxis()->SetBinLabel(3,  "EvtEmpty");          fHstatHistos->SetBinContent(3,  fNumTotEvtEmpty);
	fHstatHistos->GetXaxis()->SetBinLabel(4,  "EvtCleanEmpty");     fHstatHistos->SetBinContent(4,  fNumTotEvtCleanEmpty);
	fHstatHistos->GetXaxis()->SetBinLabel(5,  "EvtLtFem");          fHstatHistos->SetBinContent(5,  fNumTotEvtLtFem);
	fHstatHistos->GetXaxis()->SetBinLabel(6,  "EvtLtFch");          fHstatHistos->SetBinContent(6,  fNumTotEvtLtFch);
	fHstatHistos->GetXaxis()->SetBinLabel(7,  "EvtPfMETJet");       fHstatHistos->SetBinContent(7,  fNumTotEvtPfMETJet);
	fHstatHistos->GetXaxis()->SetBinLabel(8,  "EvtPfMETRij");       fHstatHistos->SetBinContent(8,  fNumTotEvtPfMETRij);
	fHstatHistos->GetXaxis()->SetBinLabel(9,  "EvtCaMETJet");       fHstatHistos->SetBinContent(9,  fNumTotEvtCaMETJet);
	fHstatHistos->GetXaxis()->SetBinLabel(10, "EvtCaMETRij");       fHstatHistos->SetBinContent(10, fNumTotEvtCaMETRij);
	fHstatHistos->GetXaxis()->SetBinLabel(11, "EvtTcMETJet");       fHstatHistos->SetBinContent(11, fNumTotEvtTcMETJet);
	fHstatHistos->GetXaxis()->SetBinLabel(12, "EvtTcMETRij");       fHstatHistos->SetBinContent(12, fNumTotEvtTcMETRij);
	fHstatHistos->GetXaxis()->SetBinLabel(13, "EvtBadHardJet");     fHstatHistos->SetBinContent(13, fNumTotEvtBadHardJet);

	fHstatHistos->GetXaxis()->SetBinLabel(21, "TotMuons");          fHstatHistos->SetBinContent(21, fNumTotMuons);
	fHstatHistos->GetXaxis()->SetBinLabel(22, "MuonGoodIso");       fHstatHistos->SetBinContent(22, fNumTotMuonGoodIso);
	fHstatHistos->GetXaxis()->SetBinLabel(23, "MuonGoodNonIso");    fHstatHistos->SetBinContent(23, fNumTotMuonGoodNonIso);
	fHstatHistos->GetXaxis()->SetBinLabel(24, "MuonBadIso");        fHstatHistos->SetBinContent(24, fNumTotMuonBadIso);
	fHstatHistos->GetXaxis()->SetBinLabel(25, "MuonBadNonIso");     fHstatHistos->SetBinContent(25, fNumTotMuonBadNonIso);
	fHstatHistos->GetXaxis()->SetBinLabel(26, "MuonDupl");          fHstatHistos->SetBinContent(26, fNumTotMuonDupl);
	fHstatHistos->GetXaxis()->SetBinLabel(27, "MuonNotPrimaryTrk"); fHstatHistos->SetBinContent(27, fNumTotMuonNotPrimaryTrk);
	fHstatHistos->GetXaxis()->SetBinLabel(28, "MuonNotClean");      fHstatHistos->SetBinContent(28, fNumTotMuonNotClean);
	fHstatHistos->GetXaxis()->SetBinLabel(29, "MuonBadDpop");       fHstatHistos->SetBinContent(29, fNumTotMuonBadDpop);
	fHstatHistos->GetXaxis()->SetBinLabel(30, "MuonBadChi2");       fHstatHistos->SetBinContent(30, fNumTotMuonBadChi2);
	fHstatHistos->GetXaxis()->SetBinLabel(31, "MuonBadNhit");       fHstatHistos->SetBinContent(31, fNumTotMuonBadNhit);

	fHstatHistos->GetXaxis()->SetBinLabel(41, "TotElectrons");      fHstatHistos->SetBinContent(41, fNumTotElectrons);
	fHstatHistos->GetXaxis()->SetBinLabel(42, "ElecGoodIso");       fHstatHistos->SetBinContent(42, fNumTotElecGoodIso);
	fHstatHistos->GetXaxis()->SetBinLabel(43, "ElecGoodNonIso");    fHstatHistos->SetBinContent(43, fNumTotElecGoodNonIso);
	fHstatHistos->GetXaxis()->SetBinLabel(44, "ElecBadIso");        fHstatHistos->SetBinContent(44, fNumTotElecBadIso);
	fHstatHistos->GetXaxis()->SetBinLabel(45, "ElecBadNonIso");     fHstatHistos->SetBinContent(45, fNumTotElecBadNonIso);
	fHstatHistos->GetXaxis()->SetBinLabel(46, "ElecDupl");          fHstatHistos->SetBinContent(46, fNumTotElecDupl);
	fHstatHistos->GetXaxis()->SetBinLabel(47, "ElecNotPrimaryTrk"); fHstatHistos->SetBinContent(47, fNumTotElecNotPrimaryTrk);
	fHstatHistos->GetXaxis()->SetBinLabel(48, "ElecNotClean");      fHstatHistos->SetBinContent(48, fNumTotElecNotClean);
	fHstatHistos->GetXaxis()->SetBinLabel(49, "ElecBadHoE");        fHstatHistos->SetBinContent(49, fNumTotElecBadHoE);
	fHstatHistos->GetXaxis()->SetBinLabel(50, "ElecBadShsh");       fHstatHistos->SetBinContent(50, fNumTotElecBadShsh);
	fHstatHistos->GetXaxis()->SetBinLabel(51, "ElecBadTmat");       fHstatHistos->SetBinContent(51, fNumTotElecBadTmat);

	fHstatHistos->GetXaxis()->SetBinLabel(61, "TotPhotons");        fHstatHistos->SetBinContent(61, fNumTotPhotons);
	fHstatHistos->GetXaxis()->SetBinLabel(62, "PhotGoodIso");       fHstatHistos->SetBinContent(62, fNumTotPhotGoodIso);
	fHstatHistos->GetXaxis()->SetBinLabel(63, "PhotGoodNonIso");    fHstatHistos->SetBinContent(63, fNumTotPhotGoodNonIso);
	fHstatHistos->GetXaxis()->SetBinLabel(64, "PhotBadIso");        fHstatHistos->SetBinContent(64, fNumTotPhotBadIso);
	fHstatHistos->GetXaxis()->SetBinLabel(65, "PhotBadNonIso");     fHstatHistos->SetBinContent(65, fNumTotPhotBadNonIso);
	fHstatHistos->GetXaxis()->SetBinLabel(66, "PhotDupl");          fHstatHistos->SetBinContent(66, fNumTotPhotDupl);
	fHstatHistos->GetXaxis()->SetBinLabel(67, "PhotNotClean");      fHstatHistos->SetBinContent(67, fNumTotPhotNotClean);
	fHstatHistos->GetXaxis()->SetBinLabel(68, "PhotBadHoE");        fHstatHistos->SetBinContent(68, fNumTotPhotBadHoE);
	fHstatHistos->GetXaxis()->SetBinLabel(69, "PhotBadShsh");       fHstatHistos->SetBinContent(69, fNumTotPhotBadShsh);

	fHstatHistos->GetXaxis()->SetBinLabel(81, "TotJets");           fHstatHistos->SetBinContent(81, fNumTotJets);
	fHstatHistos->GetXaxis()->SetBinLabel(82, "JetGood");           fHstatHistos->SetBinContent(82, fNumTotJetGood);
	fHstatHistos->GetXaxis()->SetBinLabel(83, "JetBad");            fHstatHistos->SetBinContent(83, fNumTotJetBad);
	fHstatHistos->GetXaxis()->SetBinLabel(84, "JetDuplElJet");      fHstatHistos->SetBinContent(84, fNumTotJetDuplElJet);
	fHstatHistos->GetXaxis()->SetBinLabel(85, "JetNotPrimaryTrk");  fHstatHistos->SetBinContent(85, fNumTotJetNotPrimaryTrk);
	fHstatHistos->GetXaxis()->SetBinLabel(86, "JetNotClean");       fHstatHistos->SetBinContent(86, fNumTotJetNotClean);
	fHstatHistos->GetXaxis()->SetBinLabel(87, "JetPgtE");           fHstatHistos->SetBinContent(87, fNumTotJetPgtE);
	fHstatHistos->GetXaxis()->SetBinLabel(88, "JetGtFem");          fHstatHistos->SetBinContent(88, fNumTotJetGtFem);
	fHstatHistos->GetXaxis()->SetBinLabel(89, "JetLtFem");          fHstatHistos->SetBinContent(89, fNumTotJetLtFem);
	fHstatHistos->GetXaxis()->SetBinLabel(90, "JetLtFch");          fHstatHistos->SetBinContent(90, fNumTotJetLtFch);
	fHstatHistos->GetXaxis()->SetBinLabel(91, "TotBJets");          fHstatHistos->SetBinContent(91, fNumTotBJets);

	return;
}

void TreeCleaner::PutMuon(int inew, int iold){
	// This needs to be UPDATED every time the tree content changes!
	fTR->MuGood                   [inew] = fTR->MuGood                   [iold];
	fTR->MuIsIso                  [inew] = fTR->MuIsIso                  [iold];
	fTR->MuIsGlobalMuon           [inew] = fTR->MuIsGlobalMuon           [iold];
	fTR->MuIsTrackerMuon          [inew] = fTR->MuIsTrackerMuon          [iold];
	fTR->MuPx                     [inew] = fTR->MuPx                     [iold];
	fTR->MuPy                     [inew] = fTR->MuPy                     [iold];
	fTR->MuPz                     [inew] = fTR->MuPz                     [iold];
	fTR->MuPt                     [inew] = fTR->MuPt                     [iold];
	fTR->MuPtE                    [inew] = fTR->MuPtE                    [iold];
	fTR->MuE                      [inew] = fTR->MuE                      [iold];
	fTR->MuEt                     [inew] = fTR->MuEt                     [iold];
	fTR->MuEta                    [inew] = fTR->MuEta                    [iold];
	fTR->MuPhi                    [inew] = fTR->MuPhi                    [iold];
	fTR->MuCharge                 [inew] = fTR->MuCharge                 [iold];
	fTR->MuRelIso03               [inew] = fTR->MuRelIso03               [iold];
	fTR->MuIso03SumPt             [inew] = fTR->MuIso03SumPt             [iold];
	fTR->MuIso03EmEt              [inew] = fTR->MuIso03EmEt              [iold];
	fTR->MuIso03HadEt             [inew] = fTR->MuIso03HadEt             [iold];
	fTR->MuIso03EMVetoEt          [inew] = fTR->MuIso03EMVetoEt          [iold];
	fTR->MuIso03HadVetoEt         [inew] = fTR->MuIso03HadVetoEt         [iold];
	fTR->MuIso05SumPt             [inew] = fTR->MuIso05SumPt             [iold];
	fTR->MuIso05EmEt              [inew] = fTR->MuIso05EmEt              [iold];
	fTR->MuIso05HadEt             [inew] = fTR->MuIso05HadEt             [iold];
	fTR->MuEem                    [inew] = fTR->MuEem                    [iold];
	fTR->MuEhad                   [inew] = fTR->MuEhad                   [iold];
	fTR->MuD0BS                   [inew] = fTR->MuD0BS                   [iold];
	fTR->MuD0PV                   [inew] = fTR->MuD0PV                   [iold];
	fTR->MuD0E                    [inew] = fTR->MuD0E                    [iold];
	fTR->MuDzBS                   [inew] = fTR->MuDzBS                   [iold];
	fTR->MuDzPV                   [inew] = fTR->MuDzPV                   [iold];
	fTR->MuDzE                    [inew] = fTR->MuDzE                    [iold];
	fTR->MuNChi2                  [inew] = fTR->MuNChi2                  [iold];
	fTR->MuNGlHits                [inew] = fTR->MuNGlHits                [iold];
	fTR->MuNMuHits                [inew] = fTR->MuNMuHits                [iold];
	fTR->MuNTkHits                [inew] = fTR->MuNTkHits                [iold];
	fTR->MuInnerTkNChi2           [inew] = fTR->MuInnerTkNChi2           [iold];
	fTR->MuNMatches               [inew] = fTR->MuNMatches               [iold];
	fTR->MuNChambers              [inew] = fTR->MuNChambers              [iold];
	fTR->MuCaloComp               [inew] = fTR->MuCaloComp               [iold];
	fTR->MuSegmComp               [inew] = fTR->MuSegmComp               [iold];
	fTR->MuIsGMPT                 [inew] = fTR->MuIsGMPT                 [iold];
	fTR->MuIsGMTkChiComp          [inew] = fTR->MuIsGMTkChiComp          [iold];
	fTR->MuIsGMStaChiComp         [inew] = fTR->MuIsGMStaChiComp         [iold];
	fTR->MuIsGMTkKinkTight        [inew] = fTR->MuIsGMTkKinkTight        [iold];
	fTR->MuIsAllStaMuons          [inew] = fTR->MuIsAllStaMuons          [iold];
	fTR->MuIsAllTrkMuons          [inew] = fTR->MuIsAllTrkMuons          [iold];
	fTR->MuIsTrkMuonArbitrated    [inew] = fTR->MuIsTrkMuonArbitrated    [iold];
	fTR->MuIsAllArbitrated        [inew] = fTR->MuIsAllArbitrated        [iold];
	fTR->MuIsTMLSLoose            [inew] = fTR->MuIsTMLSLoose            [iold];
	fTR->MuIsTMLSTight            [inew] = fTR->MuIsTMLSTight            [iold];
	fTR->MuIsTM2DCompLoose        [inew] = fTR->MuIsTM2DCompLoose        [iold];
	fTR->MuIsTM2DCompTight        [inew] = fTR->MuIsTM2DCompTight        [iold];
	fTR->MuIsTMOneStationLoose    [inew] = fTR->MuIsTMOneStationLoose    [iold];
	fTR->MuIsTMOneStationTight    [inew] = fTR->MuIsTMOneStationTight    [iold];
	fTR->MuIsTMLSOptLowPtLoose    [inew] = fTR->MuIsTMLSOptLowPtLoose    [iold];
	fTR->MuIsTMLSAngLoose         [inew] = fTR->MuIsTMLSAngLoose         [iold];
	fTR->MuIsTMLastStationAngTight[inew] = fTR->MuIsTMLastStationAngTight[iold];
	fTR->MuIsTMOneStationAngTight [inew] = fTR->MuIsTMOneStationAngTight [iold];
	fTR->MuIsTMOneStationAngLoose [inew] = fTR->MuIsTMOneStationAngLoose [iold];
	fTR->MuOutPosRadius           [inew] = fTR->MuOutPosRadius           [iold];
	fTR->MuOutPosX                [inew] = fTR->MuOutPosX                [iold];
	fTR->MuOutPosY                [inew] = fTR->MuOutPosY                [iold];
	fTR->MuOutPosZ                [inew] = fTR->MuOutPosZ                [iold];
	fTR->MuOutMomx                [inew] = fTR->MuOutMomx                [iold];
	fTR->MuOutMomy                [inew] = fTR->MuOutMomy                [iold];
	fTR->MuOutMomz                [inew] = fTR->MuOutMomz                [iold];
	fTR->MuOutMomPhi              [inew] = fTR->MuOutMomPhi              [iold];
	fTR->MuOutMomEta              [inew] = fTR->MuOutMomEta              [iold];
	fTR->MuOutMomTheta            [inew] = fTR->MuOutMomTheta            [iold];
	fTR->MuGenID                  [inew] = fTR->MuGenID                  [iold];
	fTR->MuGenStatus              [inew] = fTR->MuGenStatus              [iold];
	fTR->MuGenCharge              [inew] = fTR->MuGenCharge              [iold];
	fTR->MuGenPt                  [inew] = fTR->MuGenPt                  [iold];
	fTR->MuGenEta                 [inew] = fTR->MuGenEta                 [iold];
	fTR->MuGenPhi                 [inew] = fTR->MuGenPhi                 [iold];
	fTR->MuGenE                   [inew] = fTR->MuGenE                   [iold];
	fTR->MuGenMID                 [inew] = fTR->MuGenMID                 [iold];
	fTR->MuGenMStatus             [inew] = fTR->MuGenMStatus             [iold];
	fTR->MuGenMCharge             [inew] = fTR->MuGenMCharge             [iold];
	fTR->MuGenMPt                 [inew] = fTR->MuGenMPt                 [iold];
	fTR->MuGenMEta                [inew] = fTR->MuGenMEta                [iold];
	fTR->MuGenMPhi                [inew] = fTR->MuGenMPhi                [iold];
	fTR->MuGenME                  [inew] = fTR->MuGenME                  [iold];
	fTR->MuGenGMID                [inew] = fTR->MuGenGMID                [iold];
	fTR->MuGenGMStatus            [inew] = fTR->MuGenGMStatus            [iold];
	fTR->MuGenGMCharge            [inew] = fTR->MuGenGMCharge            [iold];
	fTR->MuGenGMPt                [inew] = fTR->MuGenGMPt                [iold];
	fTR->MuGenGMEta               [inew] = fTR->MuGenGMEta               [iold];
	fTR->MuGenGMPhi               [inew] = fTR->MuGenGMPhi               [iold];
	fTR->MuGenGME                 [inew] = fTR->MuGenGME                 [iold];
}

void TreeCleaner::PutElectron(int inew, int iold){
	// This needs to be UPDATED every time the tree content changes!
	fTR->ElGood                      [inew] = fTR->ElGood                     [iold];
	fTR->ElIsIso                     [inew] = fTR->ElIsIso                    [iold];
	fTR->ElChargeMisIDProb           [inew] = fTR->ElChargeMisIDProb          [iold];
	fTR->ElPx                        [inew] = fTR->ElPx                       [iold];
	fTR->ElPy                        [inew] = fTR->ElPy                       [iold];
	fTR->ElPz                        [inew] = fTR->ElPz                       [iold];
	fTR->ElPt                        [inew] = fTR->ElPt                       [iold];
	fTR->ElPtE                       [inew] = fTR->ElPtE                      [iold];
	fTR->ElE                         [inew] = fTR->ElE                        [iold];
	fTR->ElEt                        [inew] = fTR->ElEt                       [iold];
	fTR->ElEta                       [inew] = fTR->ElEta                      [iold];
	fTR->ElTheta                     [inew] = fTR->ElTheta                    [iold];
	fTR->ElSCEta                     [inew] = fTR->ElSCEta                    [iold];
	fTR->ElPhi                       [inew] = fTR->ElPhi                      [iold];
	fTR->ElD0BS                      [inew] = fTR->ElD0BS                     [iold];
	fTR->ElD0PV                      [inew] = fTR->ElD0PV                     [iold];
	fTR->ElD0E                       [inew] = fTR->ElD0E                      [iold];
	fTR->ElDzBS                      [inew] = fTR->ElDzBS                     [iold];
	fTR->ElDzPV                      [inew] = fTR->ElDzPV                     [iold];
	fTR->ElDzE                       [inew] = fTR->ElDzE                      [iold];
	// fTR->ElRelIso03                  [inew] = fTR->ElRelIso03                 [iold];
	fTR->ElRelIso04                  [inew] = fTR->ElRelIso04                 [iold];
	fTR->ElDR03TkSumPt               [inew] = fTR->ElDR03TkSumPt              [iold];
	fTR->ElDR04TkSumPt               [inew] = fTR->ElDR04TkSumPt              [iold];
	fTR->ElDR03EcalRecHitSumEt       [inew] = fTR->ElDR03EcalRecHitSumEt      [iold];
	fTR->ElDR04EcalRecHitSumEt       [inew] = fTR->ElDR04EcalRecHitSumEt      [iold];
	fTR->ElDR03HcalTowerSumEt        [inew] = fTR->ElDR03HcalTowerSumEt       [iold];
	fTR->ElDR04HcalTowerSumEt        [inew] = fTR->ElDR04HcalTowerSumEt       [iold];
	fTR->ElNChi2                     [inew] = fTR->ElNChi2                    [iold];
	fTR->ElCharge                    [inew] = fTR->ElCharge                   [iold];
	fTR->ElCInfoIsGsfCtfCons         [inew] = fTR->ElCInfoIsGsfCtfCons        [iold];
	fTR->ElCInfoIsGsfCtfScPixCons    [inew] = fTR->ElCInfoIsGsfCtfScPixCons   [iold];
	fTR->ElCInfoIsGsfScPixCons       [inew] = fTR->ElCInfoIsGsfScPixCons      [iold];
	fTR->ElScPixCharge               [inew] = fTR->ElScPixCharge              [iold];
	fTR->ElClosestCtfTrackPt         [inew] = fTR->ElClosestCtfTrackPt        [iold];
	fTR->ElClosestCtfTrackEta        [inew] = fTR->ElClosestCtfTrackEta       [iold];
	fTR->ElClosestCtfTrackPhi        [inew] = fTR->ElClosestCtfTrackPhi       [iold];
	fTR->ElClosestCtfTrackCharge     [inew] = fTR->ElClosestCtfTrackCharge    [iold];
	fTR->ElIDTight                   [inew] = fTR->ElIDTight                  [iold];
	fTR->ElIDLoose                   [inew] = fTR->ElIDLoose                  [iold];
	fTR->ElIDRobustTight             [inew] = fTR->ElIDRobustTight            [iold];
	fTR->ElIDRobustLoose             [inew] = fTR->ElIDRobustLoose            [iold];
	fTR->ElInGap                     [inew] = fTR->ElInGap                    [iold];
	fTR->ElEcalDriven                [inew] = fTR->ElEcalDriven               [iold];
	fTR->ElTrackerDriven             [inew] = fTR->ElTrackerDriven            [iold];
	fTR->ElBasicClustersSize         [inew] = fTR->ElBasicClustersSize        [iold];
	fTR->Elfbrem                     [inew] = fTR->Elfbrem                    [iold];
	fTR->ElHcalOverEcal              [inew] = fTR->ElHcalOverEcal             [iold];
	fTR->ElE1x5                      [inew] = fTR->ElE1x5                     [iold];
	fTR->ElE5x5                      [inew] = fTR->ElE5x5                     [iold];
	fTR->ElE2x5Max                   [inew] = fTR->ElE2x5Max                  [iold];
	fTR->ElSigmaIetaIeta             [inew] = fTR->ElSigmaIetaIeta            [iold];
	fTR->ElDeltaPhiSeedClusterAtCalo [inew] = fTR->ElDeltaPhiSeedClusterAtCalo[iold];
	fTR->ElDeltaEtaSeedClusterAtCalo [inew] = fTR->ElDeltaEtaSeedClusterAtCalo[iold];
	fTR->ElDeltaPhiSuperClusterAtVtx [inew] = fTR->ElDeltaPhiSuperClusterAtVtx[iold];
	fTR->ElDeltaEtaSuperClusterAtVtx [inew] = fTR->ElDeltaEtaSuperClusterAtVtx[iold];
	fTR->ElCaloEnergy                [inew] = fTR->ElCaloEnergy               [iold];
	fTR->ElTrkMomAtVtx               [inew] = fTR->ElTrkMomAtVtx              [iold];
	fTR->ElESuperClusterOverP        [inew] = fTR->ElESuperClusterOverP       [iold];
	fTR->ElNumberOfMissingInnerHits  [inew] = fTR->ElNumberOfMissingInnerHits [iold];
	fTR->ElIsInJet                   [inew] = fTR->ElIsInJet                  [iold];
	fTR->ElSharedPx                  [inew] = fTR->ElSharedPx                 [iold];
	fTR->ElSharedPy                  [inew] = fTR->ElSharedPy                 [iold];
	fTR->ElSharedPz                  [inew] = fTR->ElSharedPz                 [iold];
	fTR->ElSharedEnergy              [inew] = fTR->ElSharedEnergy             [iold];
	fTR->ElDuplicateEl               [inew] = fTR->ElDuplicateEl              [iold];
	fTR->ElConvPartnerTrkDist        [inew] = fTR->ElConvPartnerTrkDist       [iold];
	fTR->ElConvPartnerTrkDCot        [inew] = fTR->ElConvPartnerTrkDCot       [iold];
	fTR->ElConvPartnerTrkPt          [inew] = fTR->ElConvPartnerTrkPt         [iold];
	fTR->ElConvPartnerTrkEta         [inew] = fTR->ElConvPartnerTrkEta        [iold];
	fTR->ElConvPartnerTrkPhi         [inew] = fTR->ElConvPartnerTrkPhi        [iold];
	fTR->ElConvPartnerTrkCharge      [inew] = fTR->ElConvPartnerTrkCharge     [iold];
	fTR->ElScSeedSeverity            [inew] = fTR->ElScSeedSeverity           [iold];
	fTR->ElE1OverE9                  [inew] = fTR->ElE1OverE9                 [iold];
	fTR->ElS4OverS1                  [inew] = fTR->ElS4OverS1                 [iold];
	fTR->ElGenID                     [inew] = fTR->ElGenID                    [iold];
	fTR->ElGenStatus                 [inew] = fTR->ElGenStatus                [iold];
	fTR->ElGenCharge                 [inew] = fTR->ElGenCharge                [iold];
	fTR->ElGenPt                     [inew] = fTR->ElGenPt                    [iold];
	fTR->ElGenEta                    [inew] = fTR->ElGenEta                   [iold];
	fTR->ElGenPhi                    [inew] = fTR->ElGenPhi                   [iold];
	fTR->ElGenE                      [inew] = fTR->ElGenE                     [iold];
	fTR->ElGenMID                    [inew] = fTR->ElGenMID                   [iold];
	fTR->ElGenMStatus                [inew] = fTR->ElGenMStatus               [iold];
	fTR->ElGenMCharge                [inew] = fTR->ElGenMCharge               [iold];
	fTR->ElGenMPt                    [inew] = fTR->ElGenMPt                   [iold];
	fTR->ElGenMEta                   [inew] = fTR->ElGenMEta                  [iold];
	fTR->ElGenMPhi                   [inew] = fTR->ElGenMPhi                  [iold];
	fTR->ElGenME                     [inew] = fTR->ElGenME                    [iold];
	fTR->ElGenGMID                   [inew] = fTR->ElGenGMID                  [iold];
	fTR->ElGenGMStatus               [inew] = fTR->ElGenGMStatus              [iold];
	fTR->ElGenGMCharge               [inew] = fTR->ElGenGMCharge              [iold];
	fTR->ElGenGMPt                   [inew] = fTR->ElGenGMPt                  [iold];
	fTR->ElGenGMEta                  [inew] = fTR->ElGenGMEta                 [iold];
	fTR->ElGenGMPhi                  [inew] = fTR->ElGenGMPhi                 [iold];
	fTR->ElGenGME                    [inew] = fTR->ElGenGME                   [iold];
}

void TreeCleaner::PutPhoton(int inew, int iold){
	// This needs to be UPDATED every time the tree content changes!
	fTR->PhoGood          [inew] = fTR->PhoGood          [iold];
	fTR->PhoIsIso         [inew] = fTR->PhoIsIso         [iold];
	fTR->PhoPt            [inew] = fTR->PhoPt            [iold];
	fTR->PhoPx            [inew] = fTR->PhoPx            [iold];
	fTR->PhoPy            [inew] = fTR->PhoPy            [iold];
	fTR->PhoPz            [inew] = fTR->PhoPz            [iold];
	fTR->PhoEta           [inew] = fTR->PhoEta           [iold];
	fTR->PhoPhi           [inew] = fTR->PhoPhi           [iold];
	fTR->PhoEnergy        [inew] = fTR->PhoEnergy        [iold];
	fTR->PhoIso03Ecal     [inew] = fTR->PhoIso03Ecal     [iold];
	fTR->PhoIso03Hcal     [inew] = fTR->PhoIso03Hcal     [iold];
	fTR->PhoIso03TrkSolid [inew] = fTR->PhoIso03TrkSolid [iold];
	fTR->PhoIso03TrkHollow[inew] = fTR->PhoIso03TrkHollow[iold];
	fTR->PhoIso03         [inew] = fTR->PhoIso03         [iold];
	fTR->PhoCaloPositionX [inew] = fTR->PhoCaloPositionX [iold];
	fTR->PhoCaloPositionY [inew] = fTR->PhoCaloPositionY [iold];
	fTR->PhoCaloPositionZ [inew] = fTR->PhoCaloPositionZ [iold];
	fTR->PhoHoverE        [inew] = fTR->PhoHoverE        [iold];
	fTR->PhoH1overE       [inew] = fTR->PhoH1overE       [iold];
	fTR->PhoH2overE       [inew] = fTR->PhoH2overE       [iold];
	fTR->PhoSigmaIetaIeta [inew] = fTR->PhoSigmaIetaIeta [iold];
	fTR->PhoHasPixSeed    [inew] = fTR->PhoHasPixSeed    [iold];
	fTR->PhoHasConvTrks   [inew] = fTR->PhoHasConvTrks   [iold];
	fTR->PhoIsInJet       [inew] = fTR->PhoIsInJet       [iold];
	fTR->PhoIsElDupl      [inew] = fTR->PhoIsElDupl      [iold];
	fTR->PhoSharedPx      [inew] = fTR->PhoSharedPx      [iold];
	fTR->PhoSharedPy      [inew] = fTR->PhoSharedPy      [iold];
	fTR->PhoSharedPz      [inew] = fTR->PhoSharedPz      [iold];
	fTR->PhoSharedEnergy  [inew] = fTR->PhoSharedEnergy  [iold];
	fTR->PhoScSeedSeverity[inew] = fTR->PhoScSeedSeverity[iold];
	fTR->PhoE1OverE9      [inew] = fTR->PhoE1OverE9      [iold];
	fTR->PhoS4OverS1      [inew] = fTR->PhoS4OverS1      [iold];
}

void TreeCleaner::PutJet(int inew, int iold){
	// This needs to be UPDATED every time the tree content changes!
	fTR->JGood          [inew] = fTR->JGood          [iold];
	fTR->JPx            [inew] = fTR->JPx            [iold];
	fTR->JPy            [inew] = fTR->JPy            [iold];
	fTR->JPz            [inew] = fTR->JPz            [iold];
	fTR->JPt            [inew] = fTR->JPt            [iold];
	fTR->JE             [inew] = fTR->JE             [iold];
	fTR->JEt            [inew] = fTR->JEt            [iold];
	fTR->JEta           [inew] = fTR->JEta           [iold];
	fTR->JPhi           [inew] = fTR->JPhi           [iold];
	fTR->JEMfrac        [inew] = fTR->JEMfrac        [iold];
	fTR->JNConstituents [inew] = fTR->JNConstituents [iold];
	fTR->JID_HPD        [inew] = fTR->JID_HPD        [iold];
	fTR->JID_RBX        [inew] = fTR->JID_RBX        [iold];
	fTR->JID_n90Hits    [inew] = fTR->JID_n90Hits    [iold];
	fTR->JID_resEMF     [inew] = fTR->JID_resEMF     [iold];
	fTR->JID_HCALTow    [inew] = fTR->JID_HCALTow    [iold];
	fTR->JID_ECALTow    [inew] = fTR->JID_ECALTow    [iold];
	fTR->JEtaRms        [inew] = fTR->JEtaRms      [iold];
	fTR->JPhiRms        [inew] = fTR->JPhiRms      [iold];
	fTR->JbTagProbTkCntHighEff[inew]  = fTR->JbTagProbTkCntHighEff[iold];
    fTR->JbTagProbTkCntHighPur[inew]  = fTR->JbTagProbTkCntHighPur[iold];
    fTR->JbTagProbSimpSVHighEff[inew] = fTR->JbTagProbSimpSVHighEff[iold];
    fTR->JbTagProbSimpSVHighPur[inew] = fTR->JbTagProbSimpSVHighPur[iold];
	fTR->JChfrac        [inew] = fTR->JChfrac        [iold];
	fTR->JEFracHadronic [inew] = fTR->JEFracHadronic [iold];
	fTR->JMass          [inew] = fTR->JMass          [iold];
	fTR->JNAssoTracks   [inew] = fTR->JNAssoTracks   [iold];
	fTR->Jtrk1px        [inew] = fTR->Jtrk1px        [iold];
	fTR->Jtrk1py        [inew] = fTR->Jtrk1py        [iold];
	fTR->Jtrk1pz        [inew] = fTR->Jtrk1pz        [iold];
	fTR->Jtrk2px        [inew] = fTR->Jtrk2px        [iold];
	fTR->Jtrk2py        [inew] = fTR->Jtrk2py        [iold];
	fTR->Jtrk2pz        [inew] = fTR->Jtrk2pz        [iold];
	fTR->Jtrk3px        [inew] = fTR->Jtrk3px        [iold];
	fTR->Jtrk3py        [inew] = fTR->Jtrk3py        [iold];
	fTR->Jtrk3pz        [inew] = fTR->Jtrk3pz        [iold];
	fTR->JEcorr         [inew] = fTR->JEcorr         [iold];
	fTR->JeMinDR        [inew] = fTR->JeMinDR        [iold];
	fTR->JVtxx          [inew] = fTR->JVtxx          [iold];
	fTR->JVtxy          [inew] = fTR->JVtxy          [iold];
	fTR->JVtxz          [inew] = fTR->JVtxz          [iold];
	fTR->JVtxExx        [inew] = fTR->JVtxExx        [iold];
	fTR->JVtxEyx        [inew] = fTR->JVtxEyx        [iold];
	fTR->JVtxEyy        [inew] = fTR->JVtxEyy        [iold];
	fTR->JVtxEzy        [inew] = fTR->JVtxEzy        [iold];
	fTR->JVtxEzz        [inew] = fTR->JVtxEzz        [iold];
	fTR->JVtxEzx        [inew] = fTR->JVtxEzx        [iold];
	fTR->JVtxNChi2      [inew] = fTR->JVtxNChi2      [iold];
	fTR->JGenPt         [inew] = fTR->JGenPt         [iold];
	fTR->JGenEta        [inew] = fTR->JGenEta        [iold];
	fTR->JGenPhi        [inew] = fTR->JGenPhi        [iold];
	fTR->JGenE          [inew] = fTR->JGenE          [iold];
	fTR->JGenEmE        [inew] = fTR->JGenEmE        [iold];
	fTR->JGenHadE       [inew] = fTR->JGenHadE       [iold];
	fTR->JGenInvE       [inew] = fTR->JGenInvE       [iold];
}
