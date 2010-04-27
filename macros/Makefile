ROOTCFLAGS     = $(shell root-config --cflags)
ROOTLIBS       = $(shell root-config --libs)
ROOTGLIBS      = $(shell root-config --glibs)

INCLUDES       = -I./include

CXX            = g++
CXXFLAGS       = -g -fPIC -Wno-deprecated -D_GNU_SOURCE -O2 $(INCLUDES)
LD             = g++
LDFLAGS        = -g
SOFLAGS        = -shared


CXXFLAGS      += $(ROOTCFLAGS)
LIBS           = $(ROOTLIBS) 

NGLIBS         = $(ROOTGLIBS) -lMinuit -lMinuit2 -lTreePlayer
GLIBS          = $(filter-out -lNew, $(NGLIBS))

SRCS           = src/base/TreeClassBase.C src/base/TreeReader.cc src/base/TreeAnalyzerBase.cc src/base/UserAnalysisBase.cc \
                 src/UserAnalyzer.cc src/TreeAnalyzer.cc src/PhysQCAnalyzer.cc src/TreeSkimmer.cc \
                 src/UserAnalysis.cc src/DiLeptonAnalysis.cc src/TreeCleaner.cc src/MultiplicityAnalysis.cc  src/SignificanceAnalysis.cc src/PhysQCAnalysis.cc \
                 src/helper/AnaClass.cc src/helper/Davismt2.cc src/helper/LeptJetStat.cc src/helper/Hemisphere.cc src/LeptJetMultAnalysis.cc src/LeptJetMultAnalyzer.cc

OBJS           = $(patsubst %.C,%.o,$(SRCS:.cc=.o))

.SUFFIXES: .cc,.C,.hh,.h
.PHONY : clean purge all depend PhysQC

# Rules ====================================
all: RunUserAnalyzer RunTreeAnalyzer RunPhysQCAnalyzer RunTreeSkimmer RunLeptJetMultAnalyzer

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

clean:
	$(RM) $(OBJS)	
	$(RM) RunUserAnalyzer
	$(RM) RunTreeAnalyzer
	$(RM) RunPhysQCAnalyzer
	$(RM) RunTreeSkimmer
	$(RM) RunLeptJetMultAnalyzer
	

purge:
	$(RM) $(OBJS)

deps: $(SRCS)
	makedepend $(INCLUDES) $^

# DO NOT DELETE THIS LINE -- make depend needs it

