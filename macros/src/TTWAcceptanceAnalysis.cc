#include "helper/Utilities.hh"
#include "TTWAcceptanceAnalysis.hh"

using namespace std;

TTWAcceptanceAnalysis::TTWAcceptanceAnalysis(TreeReader *tr) : UserAnalysisBase(tr){
	Util::SetStyle();	
}

TTWAcceptanceAnalysis::~TTWAcceptanceAnalysis(){
}

void TTWAcceptanceAnalysis::Begin(){
	// Reset counters
	n_tot    = 0;
	n_ss     = 0;
	n_sel_mm = 0;
	n_sel_em = 0;
	n_sel_ee = 0;
	n_pp     = 0;
	n_mm     = 0;
}

void TTWAcceptanceAnalysis::Analyze(){

	// Total number of events
	n_tot++;

	// Event selection
	if(fTR->NGenLeptons < 2) return;
	
	std::vector<OrderPair> lpOrdered;
	std::vector<OrderPair> lmOrdered;
	
	for (int i = 0; i < fTR->NGenLeptons; ++i) { // loop on muons
		bool tau = ((abs(fTR->GenLeptonMID[i]) == 15) && (abs(fTR->GenLeptonGMID[i] == 24 || abs(fTR->GenLeptonGMID[i] == 6))));
		if(abs(fTR->GenLeptonMID[i]) != 6 && abs(fTR->GenLeptonMID[i]) != 24 && !tau) continue;
		if(fTR->GenLeptonID[i] ==  13) lmOrdered.push_back(make_pair(i,fTR->GenLeptonPt[i]));
		if(fTR->GenLeptonID[i] == -13) lpOrdered.push_back(make_pair(i,fTR->GenLeptonPt[i]));
	}
	for (int i = 0; i < fTR->NGenLeptons; ++i) { // loop on electrons
		bool tau = ((abs(fTR->GenLeptonMID[i]) == 15) && (abs(fTR->GenLeptonGMID[i] == 24 || abs(fTR->GenLeptonGMID[i] == 6))));
		if(abs(fTR->GenLeptonMID[i]) != 6 && abs(fTR->GenLeptonMID[i]) != 24 && !tau) continue;
		if(fTR->GenLeptonID[i] ==  11) lmOrdered.push_back(make_pair(i,fTR->GenLeptonPt[i]));
		if(fTR->GenLeptonID[i] == -11) lpOrdered.push_back(make_pair(i,fTR->GenLeptonPt[i]));
	}

	IndexByPt indexComparator; // Need this to sort collections
	std::sort(lpOrdered.begin(),lpOrdered.end(),indexComparator);
	std::sort(lmOrdered.begin(),lmOrdered.end(),indexComparator);

	// At least one same-sign pair
	if(lpOrdered.size() < 2 && lmOrdered.size() < 2) return;
	n_ss++; // Event has a same-sign pair

	std::vector< int > lp;
	std::vector< int > lm;
	if(lpOrdered.size() > 1){
		lp.push_back(lpOrdered[0].first);
		lp.push_back(lpOrdered[1].first);
	}
	if(lmOrdered.size() > 1){
		lm.push_back(lmOrdered[0].first);
		lm.push_back(lmOrdered[1].first);
	}
	
	// Select a pair
	std::vector<int> sspair;
	float sumpt1(0.), sumpt2(0.);
	int charge = 0;
	if(lp.size() > 1){
		float maxpt = std::max(fTR->GenLeptonPt[lp[0]], fTR->GenLeptonPt[lp[1]]);
		float minpt = std::min(fTR->GenLeptonPt[lp[0]], fTR->GenLeptonPt[lp[1]]);
		sumpt1 = fTR->GenLeptonPt[lp[0]] + fTR->GenLeptonPt[lp[1]];
		if(maxpt > 55. && minpt > 30.){
			sspair.push_back(lp[0]);
			sspair.push_back(lp[1]);
			charge = 1;
		}
	}
	if(lm.size() > 1){
		float maxpt = std::max(fTR->GenLeptonPt[lm[0]], fTR->GenLeptonPt[lm[1]]);
		float minpt = std::min(fTR->GenLeptonPt[lm[0]], fTR->GenLeptonPt[lm[1]]);
		sumpt2 = fTR->GenLeptonPt[lm[0]] + fTR->GenLeptonPt[lm[1]];
		if(maxpt > 55. && minpt > 30. && sumpt2 > sumpt1){
			sspair.clear();
			sspair.push_back(lm[0]);
			sspair.push_back(lm[1]);
			charge = -1;
		}
	}
		
	// Pt cuts
	if(sspair.size() < 2) return;

	// Lepton eta cuts
	if(fabs(fTR->GenLeptonEta[sspair[0]]) > 2.4) return;
	if(fabs(fTR->GenLeptonEta[sspair[1]]) > 2.4) return;

	// Flavor
	int id1 = std::abs(fTR->GenLeptonID[sspair[0]]);
	int id2 = std::abs(fTR->GenLeptonID[sspair[1]]);
	int flavor = 0;
	if(id1 == 11 && id2 == 13) flavor = 1;
	if(id1 == 13 && id2 == 11) flavor = 1;
	if(id1 == 11 && id2 == 11) flavor = 2;


	// Count jets
	int njets = 0;
	float HT = 0;
	for (int i = 0; i < fTR->NGenJets; ++i) { // loop on jets
		if(fTR->GenJetPt[i] < 20. || fabs(fTR->GenJetEta[i]) > 2.4) continue;
		if(Util::GetDeltaR(fTR->GenJetEta[i], fTR->GenLeptonEta[sspair[0]], fTR->GenJetPhi[i], fTR->GenLeptonPhi[sspair[0]]) < 0.4) continue;
		if(Util::GetDeltaR(fTR->GenJetEta[i], fTR->GenLeptonEta[sspair[1]], fTR->GenJetPhi[i], fTR->GenLeptonPhi[sspair[1]]) < 0.4) continue;
		njets++;
		HT += fTR->GenJetPt[i];
	}
	
	// Jet Selection
	if(njets < 3) return;
	if(HT < 100.) return;
	
	if(charge > 0) n_pp++;
	if(charge < 0) n_mm++;

	if(flavor == 0) n_sel_mm++;
	if(flavor == 1) n_sel_em++;
	if(flavor == 2) n_sel_ee++;
	return;
}

void TTWAcceptanceAnalysis::End(){
	long int n_tot_sel = n_sel_mm+n_sel_em+n_sel_ee;
	float ss_frac = 100*float(n_ss)/float(n_tot);
	float tot_eff = 100*float(n_tot_sel)/float(n_tot);
	float ss_eff  = 100*float(n_tot_sel)/float(n_ss);
	float av_charge = float(n_pp)/float(n_mm);
	cout << "--------------------------------------" << endl;
	cout << Form("| Tot. Events: %7d", n_tot) << endl;
	cout << Form("|   SS Events: %7d (%5.2f \%%)", n_ss, ss_frac) << endl;
	cout << Form("| Sel. Events: %7d (%5.2f \%%, %5.2f \%% of SS)", n_tot_sel, tot_eff, ss_eff) << endl;
	cout << "--------------------------------------" << endl;
	cout << "|           Selected events          |" << endl;
	cout << Form("|   mm   |   em   |   ee   |   tot   |") << endl;
	cout << Form("| %6d | %6d | %6d | %7d |", n_sel_mm, n_sel_em, n_sel_ee, n_tot_sel) << endl;
	cout << "--------------------------------------" << endl;
	cout << Form("| %5d (++), %5d (--), av = %5.2f |", n_pp, n_mm, av_charge) << endl;
	cout << "--------------------------------------" << endl;
}
