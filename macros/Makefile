#
# Macro's Makefile
#
# First run 'make deps', then: 'make [-j n] [targets]'
#
# In order to compile JES feature in UserAnalysisBase, run: 'make [-j n] "DOJES=1" [targets]'
#
ROOTCFLAGS     = $(shell root-config --cflags)
ROOTLIBS       = $(shell root-config --libs)
ROOTGLIBS      = $(shell root-config --glibs)

INCLUDES       = -I./include -I$(CMSSW_RELEASE_BASE)/src/

CXX            = g++
CXXFLAGS       = -g -fno-var-tracking -fPIC -Wno-deprecated -D_GNU_SOURCE -O2 -std=c++0x $(INCLUDES) 
LD             = g++
LDFLAGS        = -g 
SOFLAGS        = -O --no_exceptions -shared


CXXFLAGS      += $(ROOTCFLAGS)

LIBS           = $(ROOTLIBS) 


NGLIBS         = $(ROOTGLIBS) -lMinuit -lMinuit2 -lTreePlayer
GLIBS          = $(filter-out -lNew, $(NGLIBS)) 

# Juggle with patch versions...
CMSSW_BASE_VERSION = $(filter CMSSW%, $(subst _patch, , $(CMSSW_VERSION) ))
GLIBS         += -L$(CMSSW_RELEASE_BASE)/lib/$(SCRAM_ARCH) -L/swshare/cms/$(SCRAM_ARCH)/cms/cmssw/$(CMSSW_BASE_VERSION)/lib/$(SCRAM_ARCH) -lFWCoreFWLite -lFWCoreUtilities -lDataFormatsCommon -lDataFormatsFWLite -lCondFormatsJetMETObjects



SRCS           = src/base/TreeClassBase.cc src/base/TreeReader.cc src/base/TreeAnalyzerBase.cc src/base/UserAnalysisBase.cc \
                 src/UserAnalyzer.cc src/UserAnalysis.cc \
                 src/helper/PUWeight.C src/JZBAnalyzer.cc src/JZBAnalysis.cc src/QuickAnalyzer.cc src/QuickAnalysis.cc

OBJS           = $(patsubst %.C,%.o,$(SRCS:.cc=.o))

.SUFFIXES: .cc,.C,.hh,.h
.PHONY : clean purge all depend

# Rules ====================================
all: RunUserAnalyzer RunJZBAnalyzer RunQuickAnalyzer

RunUserAnalyzer: src/exe/RunUserAnalyzer.C $(OBJS)
	$(CXX) $(CXXFLAGS) -ldl $(GLIBS) $(LDFLAGS) -o $@ $^

RunJZBAnalyzer: src/exe/RunJZBAnalyzer.C $(OBJS)
	$(CXX) $(CXXFLAGS) -ldl $(GLIBS) $(LDFLAGS) -o $@ $^

RunQuickAnalyzer: src/exe/RunQuickAnalyzer.C $(OBJS)
	$(CXX) $(CXXFLAGS) -ldl $(GLIBS) $(LDFLAGS) -o $@ $^

clean:
	find src -name '*.o' -exec $(RM) -v {} ';' 
	$(RM) RunUserAnalyzer
	$(RM) RunJZBAnalyzer
	$(RM) RunQuickAnalyzer

purge:
	$(RM) $(OBJS)

deps: $(SRCS)
	makedepend $(INCLUDES) $^

# DO NOT DELETE THIS LINE -- make depend needs it
