#ifndef UserAnalysisBase_hh
#define UserAnalysisBase_hh

#include "TreeReader.hh"
#include "helper/pdgparticle.hh"
#include "helper/Utilities.hh"
#include <map>
#include <string>

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
	inline void SetOutputFile(TString file){ fOutputFile = Util::MakeOutputFile(fOutputDir + file); };

	virtual void ReadPDGTable(const char* filename = "pdgtable.txt");
	virtual int GetPDGParticle(pdgparticle&, int);
	virtual void GetHLTNames();
	virtual int GetHLTBit(string);
	virtual bool GetHLTResult(string);

	TreeReader *fTR;
	TString fOutputDir;
	TFile *fOutputFile;
	TString fTag;
	TLatex *fTlat;

	int fVerbose;
	map<int, pdgparticle> fPDGMap;	// Mapping of PDG ID names
	map<string, int> fHLTLabelMap;	// Mapping of HLT trigger bit names

	// Jet Selectors
	virtual bool IsGoodBasicJet(int);

	// Muon Selectors
	virtual bool IsGoodBasicMu(int);
	virtual bool IsTightMu(int);
	virtual bool IsLooseMu(int);
	virtual bool IsLooseNoTightMu(int);

	// Event Selectors
	virtual bool IsGoodMuEvent();
	virtual bool IsGoodElEvent();
	virtual bool SingleMuonSelection(int&);
	virtual bool DiMuonSelection(int&, int&);
	virtual bool SSDiMuonSelection(int &, int&);

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
