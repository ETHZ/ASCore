#ifndef UserAnalysisBase_hh
#define UserAnalysisBase_hh

#include "TreeReader.hh"
#include "helper/pdgparticle.hh"
#include "helper/Utilities.hh"

class UserAnalysisBase{
public:
	UserAnalysisBase(TreeReader *tr = 0);
	virtual ~UserAnalysisBase();

	virtual void Begin();
	virtual void Analyze();
	virtual void End();
	inline virtual void SetTag(TString tag){fTag = tag;};
	inline virtual void SetVerbose(int verbose){fVerbose = verbose;};

	inline void SetOutputDir(TString dir){ fOutputDir = Util::MakeOutputDir(dir); };

	virtual void ReadPDGTable(const char* filename);

	TreeReader *fTR;
	TString fOutputDir;
	TString fTag;
	TLatex *fTlat;

	int fVerbose;
	map<int, pdgparticle> fPDGMap;	// Mapping of PDG ID Names

private:

	// Object quality cuts:
	struct Cut{
		TBranch *branch;
		double upperbound;
		double lowerbound;
	};
	void ReadObjCuts(const char* = "objsel.dat");
	void ReadEvtSel(const char* = "evtsel.dat");
	bool IsGoodObj(int, std::vector<Cut>*);
	bool IsGoodEvt(std::vector<Cut>*);
	std::vector<Cut> fMuCuts;
	std::vector<Cut> fElCuts;
	std::vector<Cut> fJetCuts;
	std::vector<Cut> fEvtSelCuts;
};

#endif
