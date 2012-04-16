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
endif



SRCS           = src/base/TreeClassBase.C src/base/TreeReader.cc src/base/TreeAnalyzerBase.cc src/base/UserAnalysisBase.cc \
                 src/helper/AnaClass.cc src/helper/Davismt2.cc src/helper/FPRatios.cc src/helper/PUWeight.C src/helper/Lumi3DReWeighting_standalone.cc \
                 src/helper/FakeRatios.cc \
                 src/SSDLAnalyzer.cc src/SSDLAnalysis.cc \
                 src/SSDLPlotter.cc src/SSDLDumper.cc src/helper/MetaTreeClassBase.C

OBJS           = $(patsubst %.C,%.o,$(SRCS:.cc=.o))

.SUFFIXES: .cc,.C,.hh,.h
.PHONY : clean purge all depend

# Rules ====================================
all: RunSSDLAnalyzer RunSSDLDumper MakeSSDLPlots

RunSSDLDumper: src/exe/RunSSDLDumper.C $(OBJS)
	$(CXX) $(CXXFLAGS) -ldl $(GLIBS) $(LDFLAGS) -o $@ $^

MakeSSDLPlots: src/exe/MakeSSDLPlots.C $(OBJS)
	$(CXX) $(CXXFLAGS) -ldl $(GLIBS) $(LDFLAGS) -o $@ $^

RunSSDLAnalyzer: src/exe/RunSSDLAnalyzer.C $(OBJS)
	$(CXX) $(CXXFLAGS) -ldl $(GLIBS) $(LDFLAGS) -o $@ $^

clean:
	$(RM) $(OBJS)	
	$(RM) RunSSDLDumper
	$(RM) MakeSSDLPlots
	$(RM) RunSSDLAnalyzer

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
src/helper/FPRatios.o: ./include/helper/FPRatios.hh
src/helper/PUWeight.o: ./include/helper/PUWeight.h
src/helper/Lumi3DReWeighting_standalone.o: ./include/helper/Lumi3DReWeighting_standalone.hh
src/helper/FakeRatios.o: ./include/helper/FakeRatios.hh
src/SSDLAnalyzer.o: ./include/SSDLAnalyzer.hh
src/SSDLAnalyzer.o: ./include/base/TreeAnalyzerBase.hh
src/SSDLAnalyzer.o: ./include/base/TreeReader.hh
src/SSDLAnalyzer.o: ./include/base/TreeClassBase.h
src/SSDLAnalyzer.o: ./include/helper/Utilities.hh /usr/include/stdio.h
src/SSDLAnalyzer.o: /usr/include/features.h /usr/include/sys/cdefs.h
src/SSDLAnalyzer.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
src/SSDLAnalyzer.o: /usr/include/gnu/stubs-64.h /usr/include/bits/types.h
src/SSDLAnalyzer.o: /usr/include/bits/typesizes.h /usr/include/libio.h
src/SSDLAnalyzer.o: /usr/include/_G_config.h /usr/include/wchar.h
src/SSDLAnalyzer.o: /usr/include/bits/wchar.h /usr/include/gconv.h
src/SSDLAnalyzer.o: /usr/include/bits/stdio_lim.h
src/SSDLAnalyzer.o: /usr/include/bits/sys_errlist.h /usr/include/stdlib.h
src/SSDLAnalyzer.o: /usr/include/sys/types.h /usr/include/time.h
src/SSDLAnalyzer.o: /usr/include/endian.h /usr/include/bits/endian.h
src/SSDLAnalyzer.o: /usr/include/sys/select.h /usr/include/bits/select.h
src/SSDLAnalyzer.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
src/SSDLAnalyzer.o: /usr/include/sys/sysmacros.h
src/SSDLAnalyzer.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
src/SSDLAnalyzer.o: ./include/SSDLAnalysis.hh ./include/base/TreeReader.hh
src/SSDLAnalyzer.o: ./include/base/UserAnalysisBase.hh
src/SSDLAnalyzer.o: ./include/helper/pdgparticle.hh
src/SSDLAnalyzer.o: ./include/helper/PUWeight.h
src/SSDLAnalyzer.o: ./include/helper/Lumi3DReWeighting_standalone.hh
src/SSDLAnalyzer.o: ./include/helper/Monitor.hh
src/SSDLAnalysis.o: ./include/SSDLAnalysis.hh /usr/include/time.h
src/SSDLAnalysis.o: /usr/include/bits/types.h /usr/include/features.h
src/SSDLAnalysis.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
src/SSDLAnalysis.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h
src/SSDLAnalysis.o: /usr/include/bits/typesizes.h
src/SSDLAnalysis.o: ./include/base/TreeReader.hh
src/SSDLAnalysis.o: ./include/base/TreeClassBase.h
src/SSDLAnalysis.o: ./include/helper/Utilities.hh /usr/include/stdio.h
src/SSDLAnalysis.o: /usr/include/libio.h /usr/include/_G_config.h
src/SSDLAnalysis.o: /usr/include/wchar.h /usr/include/bits/wchar.h
src/SSDLAnalysis.o: /usr/include/gconv.h /usr/include/bits/stdio_lim.h
src/SSDLAnalysis.o: /usr/include/bits/sys_errlist.h /usr/include/stdlib.h
src/SSDLAnalysis.o: /usr/include/sys/types.h /usr/include/endian.h
src/SSDLAnalysis.o: /usr/include/bits/endian.h /usr/include/sys/select.h
src/SSDLAnalysis.o: /usr/include/bits/select.h /usr/include/bits/sigset.h
src/SSDLAnalysis.o: /usr/include/bits/time.h /usr/include/sys/sysmacros.h
src/SSDLAnalysis.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
src/SSDLAnalysis.o: ./include/base/UserAnalysisBase.hh
src/SSDLAnalysis.o: ./include/base/TreeReader.hh
src/SSDLAnalysis.o: ./include/helper/pdgparticle.hh
src/SSDLAnalysis.o: ./include/helper/PUWeight.h
src/SSDLAnalysis.o: ./include/helper/Lumi3DReWeighting_standalone.hh
src/SSDLAnalysis.o: ./include/helper/Monitor.hh
src/SSDLPlotter.o: ./include/SSDLPlotter.hh ./include/SSDLDumper.hh
src/SSDLPlotter.o: ./include/helper/AnaClass.hh /usr/include/math.h
src/SSDLPlotter.o: /usr/include/features.h /usr/include/sys/cdefs.h
src/SSDLPlotter.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
src/SSDLPlotter.o: /usr/include/gnu/stubs-64.h /usr/include/bits/huge_val.h
src/SSDLPlotter.o: /usr/include/bits/mathdef.h /usr/include/bits/mathcalls.h
src/SSDLPlotter.o: /usr/include/stdlib.h /usr/include/sys/types.h
src/SSDLPlotter.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
src/SSDLPlotter.o: /usr/include/time.h /usr/include/endian.h
src/SSDLPlotter.o: /usr/include/bits/endian.h /usr/include/sys/select.h
src/SSDLPlotter.o: /usr/include/bits/select.h /usr/include/bits/sigset.h
src/SSDLPlotter.o: /usr/include/bits/time.h /usr/include/sys/sysmacros.h
src/SSDLPlotter.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
src/SSDLPlotter.o: /usr/include/stdio.h /usr/include/libio.h
src/SSDLPlotter.o: /usr/include/_G_config.h /usr/include/wchar.h
src/SSDLPlotter.o: /usr/include/bits/wchar.h /usr/include/gconv.h
src/SSDLPlotter.o: /usr/include/bits/stdio_lim.h
src/SSDLPlotter.o: /usr/include/bits/sys_errlist.h
src/SSDLPlotter.o: ./include/helper/Utilities.hh
src/SSDLPlotter.o: ./include/helper/MetaTreeClassBase.h
src/SSDLPlotter.o: ./include/helper/Monitor.hh ./include/helper/Utilities.hh
src/SSDLPlotter.o: ./include/helper/FPRatios.hh
src/SSDLPlotter.o: ./include/helper/FakeRatios.hh
src/SSDLDumper.o: ./include/SSDLDumper.hh ./include/helper/AnaClass.hh
src/SSDLDumper.o: /usr/include/math.h /usr/include/features.h
src/SSDLDumper.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
src/SSDLDumper.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h
src/SSDLDumper.o: /usr/include/bits/huge_val.h /usr/include/bits/mathdef.h
src/SSDLDumper.o: /usr/include/bits/mathcalls.h /usr/include/stdlib.h
src/SSDLDumper.o: /usr/include/sys/types.h /usr/include/bits/types.h
src/SSDLDumper.o: /usr/include/bits/typesizes.h /usr/include/time.h
src/SSDLDumper.o: /usr/include/endian.h /usr/include/bits/endian.h
src/SSDLDumper.o: /usr/include/sys/select.h /usr/include/bits/select.h
src/SSDLDumper.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
src/SSDLDumper.o: /usr/include/sys/sysmacros.h
src/SSDLDumper.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
src/SSDLDumper.o: /usr/include/stdio.h /usr/include/libio.h
src/SSDLDumper.o: /usr/include/_G_config.h /usr/include/wchar.h
src/SSDLDumper.o: /usr/include/bits/wchar.h /usr/include/gconv.h
src/SSDLDumper.o: /usr/include/bits/stdio_lim.h
src/SSDLDumper.o: /usr/include/bits/sys_errlist.h
src/SSDLDumper.o: ./include/helper/Utilities.hh
src/SSDLDumper.o: ./include/helper/MetaTreeClassBase.h
src/SSDLDumper.o: ./include/helper/Monitor.hh ./include/helper/Utilities.hh
src/SSDLDumper.o: ./include/helper/FPRatios.hh ./include/helper/Davismt2.h
src/SSDLDumper.o: ./include/helper/FakeRatios.hh
src/helper/MetaTreeClassBase.o: ./include/helper/MetaTreeClassBase.h
