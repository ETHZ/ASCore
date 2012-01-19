run_scan(int m0){
//-----------------------------------------------------------------------------------------------------------------------------------
// pay very good attention here as to make sure you have the shared object available at the given path.
// the statistical methods can be obtained by following the instruction here: https://twiki.cern.ch/twiki/bin/view/CMS/RooStatsCl95
// make sure to have root version 5.30 or higher to compile the code.
// you can get the appropriate root version by doing `source /swshare/ROOT/thisroot.sh` or `source /swshare/ROOT/root_v5.30.00_40807_slc5_amd64//bin/thisroot.sh`
//-----------------------------------------------------------------------------------------------------------------------------------

//gROOT->ProcessLine(".L /shome/mdunser/workspace/CMSSW_4_2_8/src/DiLeptonAnalysis/NTupleProducer/macros/msugraSSDL/StatisticalTools/RooStatsRoutines/root/roostats_cl95.C");
gROOT->ProcessLine(".L StatisticalTools/RooStatsRoutines/root/roostats_cl95.C");

int   n_obs   = 4     ;
float n_exp   = 4.5   ;
float exp_err = 1.5   ;
float eff_err = 0.2   ;
float lum     = 4680. ;
float lum_err = 0.045 ;

int m0_bins(140), m0_min(220), m0_max(3000);
int m12_bins(46), m12_min(100), m12_max(1000);

//TFile * file_ = new TFile("/shome/mdunser/workspace/CMSSW_4_2_8/src/DiLeptonAnalysis/NTupleProducer/macros/msugraSSDL/res.root", "READ", "file_");
TFile * file_ = new TFile("res.root", "READ", "file_");
TH2D  * eff_  = (TH2D *) file_->Get("msugra_eff");
TH2D  * lox_  = (TH2D *) file_->Get("lo_xsec");

//TH2D *limit_     = new TH2D("msugra_limit"    , "msugra_lmit"     , m0_bins , m0_min+10, m0_max+10 , m12_bins, m12_min+10, m12_max+10);
TH2D *obs_limit_ = new TH2D("msugra_obslimit" , "msugra_obslimit" , m0_bins , m0_min-10, m0_max+10 , m12_bins, m12_min-10, m12_max+10);
TH2D *exp_limit_ = new TH2D("msugra_explimit" , "msugra_explimit" , m0_bins , m0_min-10, m0_max+10 , m12_bins, m12_min-10, m12_max+10);
TH2D *exp_hi_    = new TH2D("msugra_exphi"    , "msugra_exphi"    , m0_bins , m0_min-10, m0_max+10 , m12_bins, m12_min-10, m12_max+10);
TH2D *exp_lo_    = new TH2D("msugra_explo"    , "msugra_explo"    , m0_bins , m0_min-10, m0_max+10 , m12_bins, m12_min-10, m12_max+10);
TH2D *exp_2hi_   = new TH2D("msugra_exp2hi"   , "msugra_exp2hi"   , m0_bins , m0_min-10, m0_max+10 , m12_bins, m12_min-10, m12_max+10);
TH2D *exp_2lo_   = new TH2D("msugra_exp2lo"   , "msugra_exp2lo"   , m0_bins , m0_min-10, m0_max+10 , m12_bins, m12_min-10, m12_max+10);

bool verbose_ = true;
int m0bin = (m0-m0_min)/20 + 1;
//for (int m12 = m12_min; m12<= m12_max; m12+=20){
for (int m12bin = 1; m12bin<=m12_bins; m12bin++){
	int m12 = m12_min + (m12bin-1)*20;
	float eff = eff_->GetBinContent(m0bin, m12bin);
	if ( m0 > 1000 && m12 > 400 ) continue;
	if ( m12 > 700 ) continue;
	if ( eff == 0 ) continue;
	if ( verbose_) cout << " -----NEW POINT----------------------------------------------------------------------------------" << endl;
	if ( verbose_) cout << " ------------------------------------------------------------------------------------------------" << endl;
	if ( verbose_) cout << "m0 , m12: " << m0 << " , " << m12 << endl;
	if ( verbose_) cout << "at bins : " << m0bin << " , " << m12bin << endl;
	if ( verbose_) cout << "Efficiency: " << eff << endl;
	if ( verbose_) cout << " ------------------------------------------------------------------------------------------------" << endl;
    LimitResult limit = roostats_limit(lum, lum*lum_err, eff, eff*eff_err, n_exp, exp_err, n_obs, 0, 0, "cls", "", 0);
	if ( verbose_) cout << " ------------------------------------------------------------------------------------------------" << endl;
	if ( verbose_) cout << "Observed Limit " << limit.GetObservedLimit() << endl;
	if ( verbose_) cout << " ------------------------------------------------------------------------------------------------" << endl;
    obs_limit_ -> Fill(m0, m12, limit.GetObservedLimit());
    exp_limit_ -> Fill(m0, m12, limit.GetExpectedLimit());
    exp_hi_    -> Fill(m0, m12, limit.GetOneSigmaHighRange());
    exp_lo_    -> Fill(m0, m12, limit.GetOneSigmaLowRange());
    exp_2hi_   -> Fill(m0, m12, limit.GetTwoSigmaHighRange());
    exp_2lo_   -> Fill(m0, m12, limit.GetTwoSigmaLowRange());
}

//TFile * limits_ = new TFile(Form("/shome/mdunser/workspace/CMSSW_4_2_8/src/DiLeptonAnalysis/NTupleProducer/macros/msugraSSDL/limits_m0_%i_m12_%i_%i.root",m0, min_m12, max_m12), "CREATE", "limits_");
TFile * limits_ = new TFile(Form("limits_m0_%i_m12_%i_%i.root",m0, m12_min, m12_max), "CREATE", "limits_");
limits_->cd();
obs_limit_->Write();
exp_limit_->Write();
exp_lo_->Write();
exp_hi_->Write();
exp_2lo_->Write();
exp_2hi_->Write();

}
