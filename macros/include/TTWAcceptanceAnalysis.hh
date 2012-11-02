#ifndef TTWAcceptanceAnalysis_hh
#define TTWAcceptanceAnalysis_hh


#include <cmath>
#include <string>
#include <iostream>
#include <fstream>

#include <TH1D.h>

#include "base/TreeReader.hh"
#include "base/UserAnalysisBase.hh"

class TTWAcceptanceAnalysis : public UserAnalysisBase{
public:
	TTWAcceptanceAnalysis(TreeReader *tr = NULL);
	virtual ~TTWAcceptanceAnalysis();

	void Begin();
	void Analyze();
	void End();

	typedef std::pair<int,double> OrderPair;
	struct IndexByPt {
		const bool operator()(const OrderPair& j1, const OrderPair& j2 ) const {
			return j1.second > j2.second;
		}
	};

private:

	// Some counters
	long int n_tot;
	long int n_ss;
	long int n_sel_mm;
	long int n_sel_em;
	long int n_sel_ee;
	
	long int n_pp;
	long int n_mm;
};
#endif
