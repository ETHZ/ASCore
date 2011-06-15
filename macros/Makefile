ROOTCFLAGS     = $(shell root-config --cflags)
ROOTLIBS       = $(shell root-config --libs)
ROOTGLIBS      = $(shell root-config --glibs)

INCLUDES       = -I./include

CXX            = g++
CXXFLAGS       = -g -fPIC -Wno-deprecated -D_GNU_SOURCE -O2 $(INCLUDES) 
LD             = g++
LDFLAGS        = -g 
SOFLAGS        = -O --no_exceptions -shared


CXXFLAGS      += $(ROOTCFLAGS)
LIBS           = $(ROOTLIBS) 

NGLIBS         = $(ROOTGLIBS) -lMinuit -lMinuit2 -lTreePlayer
GLIBS          = $(filter-out -lNew, $(NGLIBS))

SRCS           = src/helper/PUWeight.C src/base/TreeClassBase.C src/base/TreeReader.cc src/base/TreeAnalyzerBase.cc src/base/UserAnalysisBase.cc \
                 src/UserAnalyzer.cc src/TreeAnalyzer.cc src/PhysQCAnalyzer.cc src/TreeSkimmer.cc src/MT2tree.cc src/MassAnalysis.cc \
                 src/UserAnalysis.cc src/DiLeptonAnalysis.cc src/TreeCleaner.cc src/MultiplicityAnalysisBase.cc \
                 src/MultiplicityAnalysis.cc  src/SignificanceAnalysis.cc src/PhysQCAnalysis.cc src/RatioAnalysis.cc \
                 src/helper/TMctLib.cc src/helper/mctlib.cc src/helper/FPRatios.cc \
                 src/helper/AnaClass.cc src/helper/Davismt2.cc src/helper/LeptJetStat.cc src/helper/Hemisphere.cc  src/LeptJetMultAnalyzer.cc \
                 src/MuonPlotter.cc src/MassPlotter.cc src/helper/MetaTreeClassBase.C \
                 src/JZBAnalyzer.cc src/JZBAnalysis.cc \
                 src/JZBPFAnalysis.cc \
                 src/SSDLAnalyzer.cc src/SSDLAnalysis.cc

# We want dictionaries only for classes that have _linkdef.h files                                                               
DICTOBS =  $(patsubst %_linkdef.hh, %.o, \
                      $(patsubst dict/%, obj/dict_%, \
                          $(wildcard dict/*_linkdef.hh) ) )

OBJS           = $(patsubst %.C,%.o,$(SRCS:.cc=.o))

OBJS += $(DICTOBS)

SHARED=shlib/libDiLeptonAnalysis.so


.SUFFIXES: .cc,.C,.hh,.h
.PHONY : clean purge all depend PhysQC

# Rules ====================================
all: RunUserAnalyzer RunTreeAnalyzer RunPhysQCAnalyzer RunTreeSkimmer RunLeptJetMultAnalyzer MakeMuonPlots RunJZBAnalyzer RunSSDLAnalyzer MakeMassPlots 
# all: RunUserAnalyzer RunTreeAnalyzer RunPhysQCAnalyzer RunTreeSkimmer RunLeptJetMultAnalyzer MakeMuonPlots RunJZBAnalyzer RunSSDLAnalyzer MakeMassPlots shared

# shared: $(SHARED)
# $(SHARED): $(OBJS)
#	@echo "Creating library $(SHARED)"
#	$(LD) $(LDFLAGS) $(SOFLAGS) $(OBJS) -o $(SHARED)
#	@echo "$(SHARED) successfully compiled!"

RunUserAnalyzer: src/exe/RunUserAnalyzer.C $(OBJS)
	$(CXX) $(CXXFLAGS) -ldl $(GLIBS) $(LDFLAGS) -o $@ $^

RunTreeAnalyzer: src/exe/RunTreeAnalyzer.C $(OBJS)
	$(CXX) $(CXXFLAGS) -ldl $(GLIBS) $(LDFLAGS) -o $@ $^

RunPhysQCAnalyzer: src/exe/RunPhysQCAnalyzer.C $(OBJS)
	$(CXX) $(CXXFLAGS) -ldl $(GLIBS) $(LDFLAGS) -o $@ $^

RunTreeSkimmer: src/exe/RunTreeSkimmer.C $(OBJS)
	$(CXX) $(CXXFLAGS) -ldl $(GLIBS) $(LDFLAGS) -o $@ $^

RunLeptJetMultAnalyzer: src/exe/RunLeptJetMultAnalyzer.C $(OBJS)
	$(CXX) $(CXXFLAGS) -ldl $(GLIBS) $(LDFLAGS) -o $@ $^	

MakeMuonPlots: src/exe/MakeMuonPlots.C $(OBJS)
	$(CXX) $(CXXFLAGS) -ldl $(GLIBS) $(LDFLAGS) -o $@ $^

MakeMassPlots: src/exe/MakeMassPlots.C $(OBJS)
	$(CXX) $(CXXFLAGS) -ldl $(GLIBS) $(LDFLAGS) -o $@ $^

RunJZBAnalyzer: src/exe/RunJZBAnalyzer.C $(OBJS)
	$(CXX) $(CXXFLAGS) -ldl $(GLIBS) $(LDFLAGS) -o $@ $^

RunSSDLAnalyzer: src/exe/RunSSDLAnalyzer.C $(OBJS)
	$(CXX) $(CXXFLAGS) -ldl $(GLIBS) $(LDFLAGS) -o $@ $^

obj/dict_%.o: include/%.hh dict/%_linkdef.hh
	@echo "Generating dictionary for $<"
	$(ROOTSYS)/bin/rootcint -f $(patsubst %.o, %.C, $@) -c -Idict $(INCLUDES) $(notdir $^)
	$(CXX) -c $(CXXFLAGS) -o $@ $(patsubst %.o, %.C, $@)

clean:
	$(RM) $(OBJS)	
	$(RM) RunUserAnalyzer
	$(RM) RunTreeAnalyzer
	$(RM) RunPhysQCAnalyzer
	$(RM) RunTreeSkimmer
	$(RM) RunLeptJetMultAnalyzer
	$(RM) MakeMuonPlots
	$(RM) MakeMassPlots
	$(RM) RunJZBAnalyzer
	$(RM) RunSSDLAnalyzer

purge:
	$(RM) $(OBJS)

deps: $(SRCS)
	makedepend $(INCLUDES) $^

# DO NOT DELETE THIS LINE -- make depend needs it

