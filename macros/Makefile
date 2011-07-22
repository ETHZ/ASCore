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

SRCS           = src/base/TreeClassBase.C src/base/TreeReader.cc src/base/TreeAnalyzerBase.cc src/base/UserAnalysisBase.cc \
                 src/UserAnalyzer.cc src/TreeSkimmer.cc src/UserAnalysis.cc \
                 src/helper/PUWeight.C src/helper/AnaClass.cc src/helper/Davismt2.cc src/helper/LeptJetStat.cc src/helper/Hemisphere.cc src/helper/MetaTreeClassBase.C

OBJS           = $(patsubst %.C,%.o,$(SRCS:.cc=.o))

.SUFFIXES: .cc,.C,.hh,.h
.PHONY : clean purge all depend

# Rules ====================================
all: RunUserAnalyzer RunTreeSkimmer 

RunUserAnalyzer: src/exe/RunUserAnalyzer.C $(OBJS)
	$(CXX) $(CXXFLAGS) -ldl $(GLIBS) $(LDFLAGS) -o $@ $^

RunTreeSkimmer: src/exe/RunTreeSkimmer.C $(OBJS)
	$(CXX) $(CXXFLAGS) -ldl $(GLIBS) $(LDFLAGS) -o $@ $^

clean:
	$(RM) $(OBJS)	
	$(RM) RunUserAnalyzer
	$(RM) RunTreeSkimmer

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
