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
CXXFLAGS       = -g -fPIC -fno-var-tracking -Wno-deprecated -D_GNU_SOURCE -O2 -std=c++0x $(INCLUDES) 
#CXXFLAGS       = -g -Wno-deprecated -D_GNU_SOURCE -O2 -std=c++0x -ftree-vectorize $(INCLUDES) 
#CXXFLAGS       = -O2 -fno-var-tracking -pedantic -ansi -pthread -pipe -Wno-vla -Werror=overflow -Wstrict-overflow -std=c++0x -msse3 -ftree-vectorize -Wno-strict-overflow -Werror=array-bounds -Werror=format-contains-nul -Werror=type-limits -fvisibility-inlines-hidden -felide-constructors -fmessage-length=0 -ftemplate-depth-300 -Wall -Wno-non-template-friend -Wno-long-long -Wreturn-type -Wunused -Wparentheses -Wno-deprecated -Werror=return-type -Werror=missing-braces -Werror=unused-value -Werror=address -Werror=format -Werror=sign-compare -Werror=write-strings -fdiagnostics-show-option -g -D_GNU_SOURCE -fPIC  $(INCLUDES)
LD             = g++
LDFLAGS        = -g 
SOFLAGS        = -O --no_exceptions -shared


CXXFLAGS      += $(ROOTCFLAGS)

LIBS           = $(ROOTLIBS) 

# Juggle with patch versions...
CMSSW_BASE_VERSION = $(filter CMSSW%, $(subst _patch, , $(CMSSW_VERSION) ))

NGLIBS         = $(ROOTGLIBS) -lMinuit -lMinuit2 -lTreePlayer
GLIBS          = $(filter-out -lNew, $(NGLIBS)) 
GLIBS         += -L$(CMSSW_RELEASE_BASE)/lib/$(SCRAM_ARCH) -L/swshare/cms/slc5_amd64_gcc462/cms/cmssw/$(CMSSW_BASE_VERSION)/lib/$(SCRAM_ARCH) -lFWCoreFWLite -lFWCoreUtilities -lDataFormatsCommon -lDataFormatsFWLite -lCondFormatsJetMETObjects



SRCS           = src/base/TreeClassBase.cc src/base/TreeReader.cc src/base/TreeAnalyzerBase.cc src/base/UserAnalysisBase.cc \
                 src/helper/PUWeight.C \
                 

OBJS           = $(patsubst %.C,%.o,$(SRCS:.cc=.o))

.SUFFIXES: .cc,.C,.hh,.h
.PHONY : clean purge all depend

# Rules ====================================
all: RunUserAnalyzer RunJZBAnalyzer #RunQuickAnalyzer

RunUserAnalyzer: src/exe/RunUserAnalyzer.C src/UserAnalyzer.cc src/UserAnalysis.cc $(OBJS)
	$(CXX) $(CXXFLAGS) -ldl $(GLIBS) $(LDFLAGS) -o $@ $^
	mv RunUserAnalyzer /scratch/$$USERNAME/RunUserAnalyzer 
	mv /scratch/$$USERNAME/RunUserAnalyzer RunUserAnalyzer

RunJZBAnalyzer: src/exe/RunJZBAnalyzer.C src/JZBAnalyzer.cc src/JZBAnalysis.cc $(OBJS)
	$(CXX) $(CXXFLAGS) -ldl $(GLIBS) $(LDFLAGS) -o $@ $^
	mv RunJZBAnalyzer /scratch/$$USERNAME/RunJZBAnalyzer
	mv /scratch/$$USERNAME/RunJZBAnalyzer RunJZBAnalyzer

#RunQuickAnalyzer: src/exe/RunQuickAnalyzer.C src/QuickAnalyzer.cc src/QuickAnalysis.cc $(OBJS)
#	$(CXX) $(CXXFLAGS) -ldl $(GLIBS) $(LDFLAGS) -o $@ $^
#	mv RunQuickAnalyzer /scratch/$$USERNAME/RunQuickAnalyzer
#	mv /scratch/$$USERNAME/RunQuickAnalyzer RunQuickAnalyzer

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


