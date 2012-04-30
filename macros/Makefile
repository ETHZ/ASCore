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
ROOFIT_INCLUDE = $(shell cd $(CMSSW_BASE); scram tool info roofitcore | grep INCLUDE= | sed 's|INCLUDE=||')
ROOFIT_LIBDIR  = $(shell cd $(CMSSW_BASE); scram tool info roofitcore | grep LIBDIR= | sed 's|LIBDIR=||')

INCLUDES       = -I./include -I$(CMSSW_RELEASE_BASE)/src/ -I$(ROOTSYS)/include  -I$(ROOFIT_INCLUDE)/

CXX            = g++
CXXFLAGS       = -g -fPIC -Wno-deprecated -D_GNU_SOURCE -O2 $(INCLUDES) 
LD             = g++
LDFLAGS        = -g 
SOFLAGS        = -O --no_exceptions -shared


CXXFLAGS      += $(ROOTCFLAGS)
ifdef DOJES
CXXFLAGS      += -D DOJES
endif

LIBS           = $(ROOTLIBS) 

NGLIBS         = $(ROOTGLIBS) -L$(ROOFIT_LIBDIR)/ -lMinuit -lMinuit2 -lTreePlayer -lRooFitCore -lRooFit
GLIBS          = $(filter-out -lNew, $(NGLIBS)) 
ifdef DOJES
GLIBS         += -L$(CMSSW_RELEASE_BASE)/lib/$(SCRAM_ARCH) -lFWCoreFWLite -lCondFormatsJetMETObjects

ifneq (,$(findstring patch,$(CMSSW_VERSION)))
	CMSSW_BASE_VERSION = $(filter CMSSW%, $(subst _patch, , $(CMSSW_VERSION) ))
	GLIBS         += -L /afs/cern.ch/cms/$(SCRAM_ARCH)/cms/cmssw/$(CMSSW_BASE_VERSION)/lib/$(SCRAM_ARCH)
endif
endif



SRCS           = src/base/TreeClassBase.C src/base/TreeReader.cc src/base/TreeAnalyzerBase.cc src/base/UserAnalysisBase.cc \
                 src/helper/AnaClass.cc src/helper/Davismt2.cc src/helper/FPRatios.cc src/helper/PUWeight.C src/helper/Lumi3DReWeighting_standalone.cc \
                 src/helper/FakeRatios.cc \
                 src/JZBAnalyzer.cc src/JZBAnalysis.cc \
                 src/helper/MetaTreeClassBase.C

OBJS           = $(patsubst %.C,%.o,$(SRCS:.cc=.o))

.SUFFIXES: .cc,.C,.hh,.h
.PHONY : clean purge all depend

# Rules ====================================
all: RunJZBAnalyzer 

RunJZBAnalyzer: src/exe/RunJZBAnalyzer.C $(OBJS) src/JZBAnalyzer.o src/JZBAnalysis.o src/JZBPFAnalysis.o
	$(CXX) $(CXXFLAGS) -ldl $(GLIBS) $(LDFLAGS) -o $@ $^

clean:
	find src -name '*.o' -exec $(RM) -v {} ';'
	$(RM) $(OBJS)	
	$(RM) RunSSDLDumper
	$(RM) MakeSSDLPlots
	$(RM) RunSSDLAnalyzer
	$(RM) RunJZBAnalyzer

purge:
	$(RM) $(OBJS)

deps: $(SRCS)
	makedepend $(INCLUDES) $^


	
# DO NOT DELETE THIS LINE -- make depend needs it
