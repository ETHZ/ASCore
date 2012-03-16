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
CXXFLAGS       = -g -fPIC -Wno-deprecated -D_GNU_SOURCE -O2 $(INCLUDES) 
LD             = g++
LDFLAGS        = -g 
SOFLAGS        = -O --no_exceptions -shared


CXXFLAGS      += $(ROOTCFLAGS)
ifdef DOJES
CXXFLAGS      += -D DOJES
endif

LIBS           = $(ROOTLIBS) 

NGLIBS         = $(ROOTGLIBS) -lMinuit -lMinuit2 -lTreePlayer
GLIBS          = $(filter-out -lNew, $(NGLIBS)) 
ifdef DOJES
GLIBS         += -L$(CMSSW_RELEASE_BASE)/lib/$(SCRAM_ARCH) -lFWCoreFWLite -lCondFormatsJetMETObjects
endif



SRCS           = src/base/TreeClassBase.C src/base/TreeReader.cc src/base/TreeAnalyzerBase.cc src/base/UserAnalysisBase.cc \
                 src/UserAnalyzer.cc src/TreeSkimmer.cc src/UserAnalysis.cc \
                 src/helper/PUWeight.C src/helper/AnaClass.cc src/helper/Davismt2.cc src/helper/LeptJetStat.cc src/helper/Hemisphere.cc src/helper/MetaTreeClassBase.C src/helper/Lumi3DReWeighting_standalone.cc src/EnergyCorrection.cc src/DiPhotonMiniTree.cc src/DiPhotonPurity.cc src/DiPhotonJetsAnalyzer.cc

OBJS           = $(patsubst %.C,%.o,$(SRCS:.cc=.o))

.SUFFIXES: .cc,.C,.hh,.h
.PHONY : clean purge all depend

# Rules ====================================
all: RunUserAnalyzer RunDiPhotonJetsAnalyzer 

RunUserAnalyzer: src/exe/RunUserAnalyzer.C $(OBJS)
	$(CXX) $(CXXFLAGS) -ldl $(GLIBS) $(LDFLAGS) -o $@ $^

RunDiPhotonJetsAnalyzer: src/exe/RunDiPhotonJetsAnalyzer.C $(OBJS)
	$(CXX) $(CXXFLAGS) -ldl $(GLIBS) $(LDFLAGS) -o $@ $^


clean:
	find src -name '*.o' -exec $(RM) -v {} ';' 
	$(RM) RunUserAnalyzer
	$(RM) RunTreeSkimmer
	$(RM) RunJZBAnalyzer
	$(RM) RunDiPhotonJetsAnalyzer

purge:
	$(RM) $(OBJS)

deps: $(SRCS)
	makedepend $(INCLUDES) $^

# DO NOT DELETE THIS LINE -- make depend needs it

src/base/TreeClassBase.o: ./include/base/TreeClassBase.h
src/base/TreeReader.o: ./include/base/TreeReader.hh
src/base/TreeReader.o: ./include/base/TreeClassBase.h
src/base/TreeAnalyzerBase.o: /usr/include/stdlib.h /usr/include/features.h
src/base/TreeAnalyzerBase.o: /usr/include/sys/cdefs.h
src/base/TreeAnalyzerBase.o: /usr/include/bits/wordsize.h
src/base/TreeAnalyzerBase.o: /usr/include/gnu/stubs.h
src/base/TreeAnalyzerBase.o: /usr/include/gnu/stubs-64.h
src/base/TreeAnalyzerBase.o: /usr/include/sys/types.h
src/base/TreeAnalyzerBase.o: /usr/include/bits/types.h
src/base/TreeAnalyzerBase.o: /usr/include/bits/typesizes.h
src/base/TreeAnalyzerBase.o: /usr/include/time.h /usr/include/endian.h
src/base/TreeAnalyzerBase.o: /usr/include/bits/endian.h
src/base/TreeAnalyzerBase.o: /usr/include/sys/select.h
src/base/TreeAnalyzerBase.o: /usr/include/bits/select.h
src/base/TreeAnalyzerBase.o: /usr/include/bits/sigset.h
src/base/TreeAnalyzerBase.o: /usr/include/bits/time.h
src/base/TreeAnalyzerBase.o: /usr/include/sys/sysmacros.h
src/base/TreeAnalyzerBase.o: /usr/include/bits/pthreadtypes.h
src/base/TreeAnalyzerBase.o: /usr/include/alloca.h
src/base/TreeAnalyzerBase.o: ./include/base/TreeAnalyzerBase.hh
src/base/TreeAnalyzerBase.o: ./include/base/TreeReader.hh
src/base/TreeAnalyzerBase.o: ./include/base/TreeClassBase.h
src/base/TreeAnalyzerBase.o: ./include/helper/Utilities.hh
src/base/TreeAnalyzerBase.o: /usr/include/stdio.h /usr/include/libio.h
src/base/TreeAnalyzerBase.o: /usr/include/_G_config.h /usr/include/wchar.h
src/base/TreeAnalyzerBase.o: /usr/include/bits/wchar.h /usr/include/gconv.h
src/base/TreeAnalyzerBase.o: /usr/include/bits/stdio_lim.h
src/base/TreeAnalyzerBase.o: /usr/include/bits/sys_errlist.h
src/base/TreeAnalyzerBase.o: /usr/include/assert.h
src/base/UserAnalysisBase.o: ./include/base/TreeReader.hh
src/base/UserAnalysisBase.o: ./include/base/TreeClassBase.h
src/base/UserAnalysisBase.o: /usr/include/stdlib.h /usr/include/features.h
src/base/UserAnalysisBase.o: /usr/include/sys/cdefs.h
src/base/UserAnalysisBase.o: /usr/include/bits/wordsize.h
src/base/UserAnalysisBase.o: /usr/include/gnu/stubs.h
src/base/UserAnalysisBase.o: /usr/include/gnu/stubs-64.h
src/base/UserAnalysisBase.o: /usr/include/sys/types.h
src/base/UserAnalysisBase.o: /usr/include/bits/types.h
src/base/UserAnalysisBase.o: /usr/include/bits/typesizes.h
src/base/UserAnalysisBase.o: /usr/include/time.h /usr/include/endian.h
src/base/UserAnalysisBase.o: /usr/include/bits/endian.h
src/base/UserAnalysisBase.o: /usr/include/sys/select.h
src/base/UserAnalysisBase.o: /usr/include/bits/select.h
src/base/UserAnalysisBase.o: /usr/include/bits/sigset.h
src/base/UserAnalysisBase.o: /usr/include/bits/time.h
src/base/UserAnalysisBase.o: /usr/include/sys/sysmacros.h
src/base/UserAnalysisBase.o: /usr/include/bits/pthreadtypes.h
src/base/UserAnalysisBase.o: /usr/include/alloca.h
src/base/UserAnalysisBase.o: ./include/helper/pdgparticle.hh
src/base/UserAnalysisBase.o: ./include/helper/Monitor.hh
src/base/UserAnalysisBase.o: ./include/base/UserAnalysisBase.hh
src/base/UserAnalysisBase.o: ./include/base/TreeReader.hh
src/base/UserAnalysisBase.o: ./include/helper/Utilities.hh
src/base/UserAnalysisBase.o: /usr/include/stdio.h /usr/include/libio.h
src/base/UserAnalysisBase.o: /usr/include/_G_config.h /usr/include/wchar.h
src/base/UserAnalysisBase.o: /usr/include/bits/wchar.h /usr/include/gconv.h
src/base/UserAnalysisBase.o: /usr/include/bits/stdio_lim.h
src/base/UserAnalysisBase.o: /usr/include/bits/sys_errlist.h
src/base/UserAnalysisBase.o: ./include/helper/PUWeight.h
src/base/UserAnalysisBase.o: ./include/helper/Lumi3DReWeighting_standalone.hh
src/UserAnalyzer.o: ./include/UserAnalyzer.hh
src/UserAnalyzer.o: ./include/base/TreeAnalyzerBase.hh
src/UserAnalyzer.o: ./include/base/TreeReader.hh
src/UserAnalyzer.o: ./include/base/TreeClassBase.h
src/UserAnalyzer.o: ./include/helper/Utilities.hh /usr/include/stdio.h
src/UserAnalyzer.o: /usr/include/features.h /usr/include/sys/cdefs.h
src/UserAnalyzer.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
src/UserAnalyzer.o: /usr/include/gnu/stubs-64.h /usr/include/bits/types.h
src/UserAnalyzer.o: /usr/include/bits/typesizes.h /usr/include/libio.h
src/UserAnalyzer.o: /usr/include/_G_config.h /usr/include/wchar.h
src/UserAnalyzer.o: /usr/include/bits/wchar.h /usr/include/gconv.h
src/UserAnalyzer.o: /usr/include/bits/stdio_lim.h
src/UserAnalyzer.o: /usr/include/bits/sys_errlist.h /usr/include/stdlib.h
src/UserAnalyzer.o: /usr/include/sys/types.h /usr/include/time.h
src/UserAnalyzer.o: /usr/include/endian.h /usr/include/bits/endian.h
src/UserAnalyzer.o: /usr/include/sys/select.h /usr/include/bits/select.h
src/UserAnalyzer.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
src/UserAnalyzer.o: /usr/include/sys/sysmacros.h
src/UserAnalyzer.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
src/UserAnalyzer.o: ./include/base/TreeReader.hh ./include/UserAnalysis.hh
src/UserAnalyzer.o: ./include/base/UserAnalysisBase.hh
src/UserAnalyzer.o: ./include/helper/pdgparticle.hh
src/UserAnalyzer.o: ./include/helper/PUWeight.h
src/UserAnalyzer.o: ./include/helper/Lumi3DReWeighting_standalone.hh
src/TreeSkimmer.o: ./include/TreeSkimmer.hh
src/TreeSkimmer.o: ./include/base/TreeAnalyzerBase.hh
src/TreeSkimmer.o: ./include/base/TreeReader.hh
src/TreeSkimmer.o: ./include/base/TreeClassBase.h
src/TreeSkimmer.o: ./include/helper/Utilities.hh /usr/include/stdio.h
src/TreeSkimmer.o: /usr/include/features.h /usr/include/sys/cdefs.h
src/TreeSkimmer.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
src/TreeSkimmer.o: /usr/include/gnu/stubs-64.h /usr/include/bits/types.h
src/TreeSkimmer.o: /usr/include/bits/typesizes.h /usr/include/libio.h
src/TreeSkimmer.o: /usr/include/_G_config.h /usr/include/wchar.h
src/TreeSkimmer.o: /usr/include/bits/wchar.h /usr/include/gconv.h
src/TreeSkimmer.o: /usr/include/bits/stdio_lim.h
src/TreeSkimmer.o: /usr/include/bits/sys_errlist.h /usr/include/stdlib.h
src/TreeSkimmer.o: /usr/include/sys/types.h /usr/include/time.h
src/TreeSkimmer.o: /usr/include/endian.h /usr/include/bits/endian.h
src/TreeSkimmer.o: /usr/include/sys/select.h /usr/include/bits/select.h
src/TreeSkimmer.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
src/TreeSkimmer.o: /usr/include/sys/sysmacros.h
src/TreeSkimmer.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
src/TreeSkimmer.o: ./include/base/TreeReader.hh
src/UserAnalysis.o: ./include/helper/Utilities.hh /usr/include/stdio.h
src/UserAnalysis.o: /usr/include/features.h /usr/include/sys/cdefs.h
src/UserAnalysis.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
src/UserAnalysis.o: /usr/include/gnu/stubs-64.h /usr/include/bits/types.h
src/UserAnalysis.o: /usr/include/bits/typesizes.h /usr/include/libio.h
src/UserAnalysis.o: /usr/include/_G_config.h /usr/include/wchar.h
src/UserAnalysis.o: /usr/include/bits/wchar.h /usr/include/gconv.h
src/UserAnalysis.o: /usr/include/bits/stdio_lim.h
src/UserAnalysis.o: /usr/include/bits/sys_errlist.h /usr/include/stdlib.h
src/UserAnalysis.o: /usr/include/sys/types.h /usr/include/time.h
src/UserAnalysis.o: /usr/include/endian.h /usr/include/bits/endian.h
src/UserAnalysis.o: /usr/include/sys/select.h /usr/include/bits/select.h
src/UserAnalysis.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
src/UserAnalysis.o: /usr/include/sys/sysmacros.h
src/UserAnalysis.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
src/UserAnalysis.o: ./include/UserAnalysis.hh ./include/base/TreeReader.hh
src/UserAnalysis.o: ./include/base/TreeClassBase.h
src/UserAnalysis.o: ./include/base/UserAnalysisBase.hh
src/UserAnalysis.o: ./include/base/TreeReader.hh
src/UserAnalysis.o: ./include/helper/pdgparticle.hh
src/UserAnalysis.o: ./include/helper/PUWeight.h
src/UserAnalysis.o: ./include/helper/Lumi3DReWeighting_standalone.hh
src/helper/PUWeight.o: ./include/helper/PUWeight.h
src/helper/AnaClass.o: ./include/helper/AnaClass.hh /usr/include/math.h
src/helper/AnaClass.o: /usr/include/features.h /usr/include/sys/cdefs.h
src/helper/AnaClass.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
src/helper/AnaClass.o: /usr/include/gnu/stubs-64.h
src/helper/AnaClass.o: /usr/include/bits/huge_val.h
src/helper/AnaClass.o: /usr/include/bits/mathdef.h
src/helper/AnaClass.o: /usr/include/bits/mathcalls.h /usr/include/stdlib.h
src/helper/AnaClass.o: /usr/include/sys/types.h /usr/include/bits/types.h
src/helper/AnaClass.o: /usr/include/bits/typesizes.h /usr/include/time.h
src/helper/AnaClass.o: /usr/include/endian.h /usr/include/bits/endian.h
src/helper/AnaClass.o: /usr/include/sys/select.h /usr/include/bits/select.h
src/helper/AnaClass.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
src/helper/AnaClass.o: /usr/include/sys/sysmacros.h
src/helper/AnaClass.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
src/helper/AnaClass.o: /usr/include/stdio.h /usr/include/libio.h
src/helper/AnaClass.o: /usr/include/_G_config.h /usr/include/wchar.h
src/helper/AnaClass.o: /usr/include/bits/wchar.h /usr/include/gconv.h
src/helper/AnaClass.o: /usr/include/bits/stdio_lim.h
src/helper/AnaClass.o: /usr/include/bits/sys_errlist.h
src/helper/AnaClass.o: ./include/helper/Utilities.hh
src/helper/AnaClass.o: ./include/helper/MetaTreeClassBase.h
src/helper/AnaClass.o: ./include/helper/Utilities.hh
src/helper/Davismt2.o: ./include/helper/Davismt2.h /usr/include/math.h
src/helper/Davismt2.o: /usr/include/features.h /usr/include/sys/cdefs.h
src/helper/Davismt2.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
src/helper/Davismt2.o: /usr/include/gnu/stubs-64.h
src/helper/Davismt2.o: /usr/include/bits/huge_val.h
src/helper/Davismt2.o: /usr/include/bits/mathdef.h
src/helper/Davismt2.o: /usr/include/bits/mathcalls.h
src/helper/LeptJetStat.o: ./include/helper/LeptJetStat.h
src/helper/Hemisphere.o: ./include/helper/Hemisphere.hh
src/helper/Hemisphere.o: ./include/helper/Utilities.hh /usr/include/stdio.h
src/helper/Hemisphere.o: /usr/include/features.h /usr/include/sys/cdefs.h
src/helper/Hemisphere.o: /usr/include/bits/wordsize.h
src/helper/Hemisphere.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h
src/helper/Hemisphere.o: /usr/include/bits/types.h
src/helper/Hemisphere.o: /usr/include/bits/typesizes.h /usr/include/libio.h
src/helper/Hemisphere.o: /usr/include/_G_config.h /usr/include/wchar.h
src/helper/Hemisphere.o: /usr/include/bits/wchar.h /usr/include/gconv.h
src/helper/Hemisphere.o: /usr/include/bits/stdio_lim.h
src/helper/Hemisphere.o: /usr/include/bits/sys_errlist.h
src/helper/Hemisphere.o: /usr/include/stdlib.h /usr/include/sys/types.h
src/helper/Hemisphere.o: /usr/include/time.h /usr/include/endian.h
src/helper/Hemisphere.o: /usr/include/bits/endian.h /usr/include/sys/select.h
src/helper/Hemisphere.o: /usr/include/bits/select.h
src/helper/Hemisphere.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
src/helper/Hemisphere.o: /usr/include/sys/sysmacros.h
src/helper/Hemisphere.o: /usr/include/bits/pthreadtypes.h
src/helper/Hemisphere.o: /usr/include/alloca.h
src/helper/MetaTreeClassBase.o: ./include/helper/MetaTreeClassBase.h
src/helper/Lumi3DReWeighting_standalone.o: ./include/helper/Lumi3DReWeighting_standalone.hh
src/EnergyCorrection.o: ./include/EnergyCorrection.hh
src/EnergyCorrection.o: ./include/base/TreeReader.hh
src/EnergyCorrection.o: ./include/base/TreeClassBase.h /usr/include/assert.h
src/EnergyCorrection.o: /usr/include/features.h /usr/include/sys/cdefs.h
src/EnergyCorrection.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
src/EnergyCorrection.o: /usr/include/gnu/stubs-64.h
src/DiPhotonMiniTree.o: ./include/helper/Utilities.hh /usr/include/stdio.h
src/DiPhotonMiniTree.o: /usr/include/features.h /usr/include/sys/cdefs.h
src/DiPhotonMiniTree.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
src/DiPhotonMiniTree.o: /usr/include/gnu/stubs-64.h /usr/include/bits/types.h
src/DiPhotonMiniTree.o: /usr/include/bits/typesizes.h /usr/include/libio.h
src/DiPhotonMiniTree.o: /usr/include/_G_config.h /usr/include/wchar.h
src/DiPhotonMiniTree.o: /usr/include/bits/wchar.h /usr/include/gconv.h
src/DiPhotonMiniTree.o: /usr/include/bits/stdio_lim.h
src/DiPhotonMiniTree.o: /usr/include/bits/sys_errlist.h /usr/include/stdlib.h
src/DiPhotonMiniTree.o: /usr/include/sys/types.h /usr/include/time.h
src/DiPhotonMiniTree.o: /usr/include/endian.h /usr/include/bits/endian.h
src/DiPhotonMiniTree.o: /usr/include/sys/select.h /usr/include/bits/select.h
src/DiPhotonMiniTree.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
src/DiPhotonMiniTree.o: /usr/include/sys/sysmacros.h
src/DiPhotonMiniTree.o: /usr/include/bits/pthreadtypes.h
src/DiPhotonMiniTree.o: /usr/include/alloca.h ./include/DiPhotonMiniTree.hh
src/DiPhotonMiniTree.o: ./include/base/TreeReader.hh
src/DiPhotonMiniTree.o: ./include/base/TreeClassBase.h
src/DiPhotonMiniTree.o: ./include/base/UserAnalysisBase.hh
src/DiPhotonMiniTree.o: ./include/base/TreeReader.hh
src/DiPhotonMiniTree.o: ./include/helper/pdgparticle.hh
src/DiPhotonMiniTree.o: ./include/helper/PUWeight.h
src/DiPhotonMiniTree.o: ./include/helper/Lumi3DReWeighting_standalone.hh
src/DiPhotonMiniTree.o: ./include/EnergyCorrection.hh
src/DiPhotonMiniTree.o: ./include/DiPhotonPurity.hh /usr/include/assert.h
src/DiPhotonMiniTree.o: src/DiPhotonTriggerSelection.cc
src/DiPhotonPurity.o: ./include/helper/Utilities.hh /usr/include/stdio.h
src/DiPhotonPurity.o: /usr/include/features.h /usr/include/sys/cdefs.h
src/DiPhotonPurity.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
src/DiPhotonPurity.o: /usr/include/gnu/stubs-64.h /usr/include/bits/types.h
src/DiPhotonPurity.o: /usr/include/bits/typesizes.h /usr/include/libio.h
src/DiPhotonPurity.o: /usr/include/_G_config.h /usr/include/wchar.h
src/DiPhotonPurity.o: /usr/include/bits/wchar.h /usr/include/gconv.h
src/DiPhotonPurity.o: /usr/include/bits/stdio_lim.h
src/DiPhotonPurity.o: /usr/include/bits/sys_errlist.h /usr/include/stdlib.h
src/DiPhotonPurity.o: /usr/include/sys/types.h /usr/include/time.h
src/DiPhotonPurity.o: /usr/include/endian.h /usr/include/bits/endian.h
src/DiPhotonPurity.o: /usr/include/sys/select.h /usr/include/bits/select.h
src/DiPhotonPurity.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
src/DiPhotonPurity.o: /usr/include/sys/sysmacros.h
src/DiPhotonPurity.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
src/DiPhotonPurity.o: ./include/DiPhotonPurity.hh
src/DiPhotonPurity.o: ./include/base/TreeReader.hh
src/DiPhotonPurity.o: ./include/base/TreeClassBase.h
src/DiPhotonPurity.o: ./include/base/UserAnalysisBase.hh
src/DiPhotonPurity.o: ./include/base/TreeReader.hh
src/DiPhotonPurity.o: ./include/helper/pdgparticle.hh
src/DiPhotonPurity.o: ./include/helper/PUWeight.h
src/DiPhotonPurity.o: ./include/helper/Lumi3DReWeighting_standalone.hh
src/DiPhotonPurity.o: ./include/EnergyCorrection.hh /usr/include/assert.h
src/DiPhotonPurity.o: src/DiPhotonTriggerSelection.cc
src/DiPhotonJetsAnalyzer.o: ./include/DiPhotonJetsAnalyzer.hh
src/DiPhotonJetsAnalyzer.o: ./include/base/TreeAnalyzerBase.hh
src/DiPhotonJetsAnalyzer.o: ./include/base/TreeReader.hh
src/DiPhotonJetsAnalyzer.o: ./include/base/TreeClassBase.h
src/DiPhotonJetsAnalyzer.o: ./include/helper/Utilities.hh
src/DiPhotonJetsAnalyzer.o: /usr/include/stdio.h /usr/include/features.h
src/DiPhotonJetsAnalyzer.o: /usr/include/sys/cdefs.h
src/DiPhotonJetsAnalyzer.o: /usr/include/bits/wordsize.h
src/DiPhotonJetsAnalyzer.o: /usr/include/gnu/stubs.h
src/DiPhotonJetsAnalyzer.o: /usr/include/gnu/stubs-64.h
src/DiPhotonJetsAnalyzer.o: /usr/include/bits/types.h
src/DiPhotonJetsAnalyzer.o: /usr/include/bits/typesizes.h
src/DiPhotonJetsAnalyzer.o: /usr/include/libio.h /usr/include/_G_config.h
src/DiPhotonJetsAnalyzer.o: /usr/include/wchar.h /usr/include/bits/wchar.h
src/DiPhotonJetsAnalyzer.o: /usr/include/gconv.h
src/DiPhotonJetsAnalyzer.o: /usr/include/bits/stdio_lim.h
src/DiPhotonJetsAnalyzer.o: /usr/include/bits/sys_errlist.h
src/DiPhotonJetsAnalyzer.o: /usr/include/stdlib.h /usr/include/sys/types.h
src/DiPhotonJetsAnalyzer.o: /usr/include/time.h /usr/include/endian.h
src/DiPhotonJetsAnalyzer.o: /usr/include/bits/endian.h
src/DiPhotonJetsAnalyzer.o: /usr/include/sys/select.h
src/DiPhotonJetsAnalyzer.o: /usr/include/bits/select.h
src/DiPhotonJetsAnalyzer.o: /usr/include/bits/sigset.h
src/DiPhotonJetsAnalyzer.o: /usr/include/bits/time.h
src/DiPhotonJetsAnalyzer.o: /usr/include/sys/sysmacros.h
src/DiPhotonJetsAnalyzer.o: /usr/include/bits/pthreadtypes.h
src/DiPhotonJetsAnalyzer.o: /usr/include/alloca.h
src/DiPhotonJetsAnalyzer.o: ./include/base/TreeReader.hh
src/DiPhotonJetsAnalyzer.o: ./include/DiPhotonPurity.hh
src/DiPhotonJetsAnalyzer.o: ./include/base/UserAnalysisBase.hh
src/DiPhotonJetsAnalyzer.o: ./include/helper/pdgparticle.hh
src/DiPhotonJetsAnalyzer.o: ./include/helper/PUWeight.h
src/DiPhotonJetsAnalyzer.o: ./include/helper/Lumi3DReWeighting_standalone.hh
src/DiPhotonJetsAnalyzer.o: ./include/EnergyCorrection.hh
src/DiPhotonJetsAnalyzer.o: ./include/DiPhotonMiniTree.hh
