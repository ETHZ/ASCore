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
                 src/UserAnalyzer.cc src/TreeAnalyzer.cc src/PhysQCAnalyzer.cc src/TreeSkimmer.cc src/MassAnalysis.cc \
                 src/UserAnalysis.cc src/DiLeptonAnalysis.cc src/TreeCleaner.cc src/MultiplicityAnalysisBase.cc \
                 src/MultiplicityAnalysis.cc  src/SignificanceAnalysis.cc src/PhysQCAnalysis.cc src/RatioAnalysis.cc \
                 src/helper/AnaClass.cc src/helper/Davismt2.cc src/helper/LeptJetStat.cc src/helper/Hemisphere.cc  src/LeptJetMultAnalyzer.cc \
                 src/MuonFakeRatioFiller.cc src/MuonFakeAnalyzer.cc src/MuonPlotter.cc \
                 src/helper/FPRatios.cc \

OBJS           = $(patsubst %.C,%.o,$(SRCS:.cc=.o))

.SUFFIXES: .cc,.C,.hh,.h
.PHONY : clean purge all depend PhysQC

# Rules ====================================
all: RunUserAnalyzer RunTreeAnalyzer RunPhysQCAnalyzer RunTreeSkimmer RunLeptJetMultAnalyzer RunMuonFakeAnalyzer MakeMuonPlots

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

RunMuonFakeAnalyzer: src/exe/RunMuonFakeAnalyzer.C $(OBJS)
	$(CXX) $(CXXFLAGS) -ldl $(GLIBS) $(LDFLAGS) -o $@ $^


clean:
	$(RM) $(OBJS)	
	$(RM) RunUserAnalyzer
	$(RM) RunTreeAnalyzer
	$(RM) RunPhysQCAnalyzer
	$(RM) RunTreeSkimmer
	$(RM) RunLeptJetMultAnalyzer
	$(RM) RunMuonFakeAnalyzer
	$(RM) MakeMuonPlots

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
src/base/UserAnalysisBase.o: ./include/base/UserAnalysisBase.hh
src/base/UserAnalysisBase.o: ./include/base/TreeReader.hh
src/base/UserAnalysisBase.o: ./include/helper/Utilities.hh
src/base/UserAnalysisBase.o: /usr/include/stdio.h /usr/include/libio.h
src/base/UserAnalysisBase.o: /usr/include/_G_config.h /usr/include/wchar.h
src/base/UserAnalysisBase.o: /usr/include/bits/wchar.h /usr/include/gconv.h
src/base/UserAnalysisBase.o: /usr/include/bits/stdio_lim.h
src/base/UserAnalysisBase.o: /usr/include/bits/sys_errlist.h
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
src/TreeAnalyzer.o: ./include/TreeAnalyzer.hh
src/TreeAnalyzer.o: ./include/base/TreeAnalyzerBase.hh
src/TreeAnalyzer.o: ./include/base/TreeReader.hh
src/TreeAnalyzer.o: ./include/base/TreeClassBase.h
src/TreeAnalyzer.o: ./include/helper/Utilities.hh /usr/include/stdio.h
src/TreeAnalyzer.o: /usr/include/features.h /usr/include/sys/cdefs.h
src/TreeAnalyzer.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
src/TreeAnalyzer.o: /usr/include/gnu/stubs-64.h /usr/include/bits/types.h
src/TreeAnalyzer.o: /usr/include/bits/typesizes.h /usr/include/libio.h
src/TreeAnalyzer.o: /usr/include/_G_config.h /usr/include/wchar.h
src/TreeAnalyzer.o: /usr/include/bits/wchar.h /usr/include/gconv.h
src/TreeAnalyzer.o: /usr/include/bits/stdio_lim.h
src/TreeAnalyzer.o: /usr/include/bits/sys_errlist.h /usr/include/stdlib.h
src/TreeAnalyzer.o: /usr/include/sys/types.h /usr/include/time.h
src/TreeAnalyzer.o: /usr/include/endian.h /usr/include/bits/endian.h
src/TreeAnalyzer.o: /usr/include/sys/select.h /usr/include/bits/select.h
src/TreeAnalyzer.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
src/TreeAnalyzer.o: /usr/include/sys/sysmacros.h
src/TreeAnalyzer.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
src/TreeAnalyzer.o: ./include/base/TreeReader.hh ./include/TreeCleaner.hh
src/TreeAnalyzer.o: ./include/base/UserAnalysisBase.hh
src/TreeAnalyzer.o: ./include/helper/pdgparticle.hh
src/TreeAnalyzer.o: ./include/helper/Davismt2.h /usr/include/math.h
src/TreeAnalyzer.o: /usr/include/bits/huge_val.h /usr/include/bits/mathdef.h
src/TreeAnalyzer.o: /usr/include/bits/mathcalls.h
src/TreeAnalyzer.o: ./include/DiLeptonAnalysis.hh
src/TreeAnalyzer.o: ./include/MultiplicityAnalysis.hh
src/TreeAnalyzer.o: ./include/MultiplicityAnalysisBase.hh
src/TreeAnalyzer.o: ./include/helper/LeptJetStat.h
src/TreeAnalyzer.o: ./include/SignificanceAnalysis.hh
src/PhysQCAnalyzer.o: ./include/PhysQCAnalyzer.hh
src/PhysQCAnalyzer.o: ./include/base/TreeAnalyzerBase.hh
src/PhysQCAnalyzer.o: ./include/base/TreeReader.hh
src/PhysQCAnalyzer.o: ./include/base/TreeClassBase.h
src/PhysQCAnalyzer.o: ./include/helper/Utilities.hh /usr/include/stdio.h
src/PhysQCAnalyzer.o: /usr/include/features.h /usr/include/sys/cdefs.h
src/PhysQCAnalyzer.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
src/PhysQCAnalyzer.o: /usr/include/gnu/stubs-64.h /usr/include/bits/types.h
src/PhysQCAnalyzer.o: /usr/include/bits/typesizes.h /usr/include/libio.h
src/PhysQCAnalyzer.o: /usr/include/_G_config.h /usr/include/wchar.h
src/PhysQCAnalyzer.o: /usr/include/bits/wchar.h /usr/include/gconv.h
src/PhysQCAnalyzer.o: /usr/include/bits/stdio_lim.h
src/PhysQCAnalyzer.o: /usr/include/bits/sys_errlist.h /usr/include/stdlib.h
src/PhysQCAnalyzer.o: /usr/include/sys/types.h /usr/include/time.h
src/PhysQCAnalyzer.o: /usr/include/endian.h /usr/include/bits/endian.h
src/PhysQCAnalyzer.o: /usr/include/sys/select.h /usr/include/bits/select.h
src/PhysQCAnalyzer.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
src/PhysQCAnalyzer.o: /usr/include/sys/sysmacros.h
src/PhysQCAnalyzer.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
src/PhysQCAnalyzer.o: ./include/base/TreeReader.hh ./include/TreeCleaner.hh
src/PhysQCAnalyzer.o: ./include/base/UserAnalysisBase.hh
src/PhysQCAnalyzer.o: ./include/helper/pdgparticle.hh
src/PhysQCAnalyzer.o: ./include/helper/Davismt2.h /usr/include/math.h
src/PhysQCAnalyzer.o: /usr/include/bits/huge_val.h
src/PhysQCAnalyzer.o: /usr/include/bits/mathdef.h
src/PhysQCAnalyzer.o: /usr/include/bits/mathcalls.h
src/PhysQCAnalyzer.o: ./include/PhysQCAnalysis.hh
src/PhysQCAnalyzer.o: ./include/helper/AnaClass.hh
src/PhysQCAnalyzer.o: ./include/helper/Utilities.hh
src/PhysQCAnalyzer.o: ./include/MultiplicityAnalysis.hh
src/PhysQCAnalyzer.o: ./include/MultiplicityAnalysisBase.hh
src/PhysQCAnalyzer.o: ./include/helper/LeptJetStat.h
src/PhysQCAnalyzer.o: ./include/DiLeptonAnalysis.hh
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
src/MassAnalysis.o: ./include/helper/Utilities.hh /usr/include/stdio.h
src/MassAnalysis.o: /usr/include/features.h /usr/include/sys/cdefs.h
src/MassAnalysis.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
src/MassAnalysis.o: /usr/include/gnu/stubs-64.h /usr/include/bits/types.h
src/MassAnalysis.o: /usr/include/bits/typesizes.h /usr/include/libio.h
src/MassAnalysis.o: /usr/include/_G_config.h /usr/include/wchar.h
src/MassAnalysis.o: /usr/include/bits/wchar.h /usr/include/gconv.h
src/MassAnalysis.o: /usr/include/bits/stdio_lim.h
src/MassAnalysis.o: /usr/include/bits/sys_errlist.h /usr/include/stdlib.h
src/MassAnalysis.o: /usr/include/sys/types.h /usr/include/time.h
src/MassAnalysis.o: /usr/include/endian.h /usr/include/bits/endian.h
src/MassAnalysis.o: /usr/include/sys/select.h /usr/include/bits/select.h
src/MassAnalysis.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
src/MassAnalysis.o: /usr/include/sys/sysmacros.h
src/MassAnalysis.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
src/MassAnalysis.o: ./include/MassAnalysis.hh ./include/base/TreeReader.hh
src/MassAnalysis.o: ./include/base/TreeClassBase.h
src/MassAnalysis.o: ./include/MultiplicityAnalysisBase.hh
src/MassAnalysis.o: ./include/base/UserAnalysisBase.hh
src/MassAnalysis.o: ./include/base/TreeReader.hh
src/MassAnalysis.o: ./include/helper/pdgparticle.hh
src/MassAnalysis.o: ./include/helper/Davismt2.h /usr/include/math.h
src/MassAnalysis.o: /usr/include/bits/huge_val.h /usr/include/bits/mathdef.h
src/MassAnalysis.o: /usr/include/bits/mathcalls.h
src/MassAnalysis.o: ./include/helper/Hemisphere.hh
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
src/DiLeptonAnalysis.o: ./include/DiLeptonAnalysis.hh
src/DiLeptonAnalysis.o: ./include/base/TreeReader.hh
src/DiLeptonAnalysis.o: ./include/base/TreeClassBase.h
src/DiLeptonAnalysis.o: ./include/base/UserAnalysisBase.hh
src/DiLeptonAnalysis.o: ./include/base/TreeReader.hh
src/DiLeptonAnalysis.o: ./include/helper/pdgparticle.hh /usr/include/stdlib.h
src/DiLeptonAnalysis.o: /usr/include/features.h /usr/include/sys/cdefs.h
src/DiLeptonAnalysis.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
src/DiLeptonAnalysis.o: /usr/include/gnu/stubs-64.h /usr/include/sys/types.h
src/DiLeptonAnalysis.o: /usr/include/bits/types.h
src/DiLeptonAnalysis.o: /usr/include/bits/typesizes.h /usr/include/time.h
src/DiLeptonAnalysis.o: /usr/include/endian.h /usr/include/bits/endian.h
src/DiLeptonAnalysis.o: /usr/include/sys/select.h /usr/include/bits/select.h
src/DiLeptonAnalysis.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
src/DiLeptonAnalysis.o: /usr/include/sys/sysmacros.h
src/DiLeptonAnalysis.o: /usr/include/bits/pthreadtypes.h
src/DiLeptonAnalysis.o: /usr/include/alloca.h ./include/helper/Utilities.hh
src/DiLeptonAnalysis.o: /usr/include/stdio.h /usr/include/libio.h
src/DiLeptonAnalysis.o: /usr/include/_G_config.h /usr/include/wchar.h
src/DiLeptonAnalysis.o: /usr/include/bits/wchar.h /usr/include/gconv.h
src/DiLeptonAnalysis.o: /usr/include/bits/stdio_lim.h
src/DiLeptonAnalysis.o: /usr/include/bits/sys_errlist.h
src/DiLeptonAnalysis.o: ./include/helper/Davismt2.h /usr/include/math.h
src/DiLeptonAnalysis.o: /usr/include/bits/huge_val.h
src/DiLeptonAnalysis.o: /usr/include/bits/mathdef.h
src/DiLeptonAnalysis.o: /usr/include/bits/mathcalls.h
src/TreeCleaner.o: ./include/TreeCleaner.hh ./include/base/TreeReader.hh
src/TreeCleaner.o: ./include/base/TreeClassBase.h
src/TreeCleaner.o: ./include/base/UserAnalysisBase.hh
src/TreeCleaner.o: ./include/base/TreeReader.hh
src/TreeCleaner.o: ./include/helper/pdgparticle.hh /usr/include/stdlib.h
src/TreeCleaner.o: /usr/include/features.h /usr/include/sys/cdefs.h
src/TreeCleaner.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
src/TreeCleaner.o: /usr/include/gnu/stubs-64.h /usr/include/sys/types.h
src/TreeCleaner.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
src/TreeCleaner.o: /usr/include/time.h /usr/include/endian.h
src/TreeCleaner.o: /usr/include/bits/endian.h /usr/include/sys/select.h
src/TreeCleaner.o: /usr/include/bits/select.h /usr/include/bits/sigset.h
src/TreeCleaner.o: /usr/include/bits/time.h /usr/include/sys/sysmacros.h
src/TreeCleaner.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
src/TreeCleaner.o: ./include/helper/Utilities.hh /usr/include/stdio.h
src/TreeCleaner.o: /usr/include/libio.h /usr/include/_G_config.h
src/TreeCleaner.o: /usr/include/wchar.h /usr/include/bits/wchar.h
src/TreeCleaner.o: /usr/include/gconv.h /usr/include/bits/stdio_lim.h
src/TreeCleaner.o: /usr/include/bits/sys_errlist.h
src/TreeCleaner.o: ./include/helper/Davismt2.h /usr/include/math.h
src/TreeCleaner.o: /usr/include/bits/huge_val.h /usr/include/bits/mathdef.h
src/TreeCleaner.o: /usr/include/bits/mathcalls.h
src/MultiplicityAnalysisBase.o: ./include/helper/Utilities.hh
src/MultiplicityAnalysisBase.o: /usr/include/stdio.h /usr/include/features.h
src/MultiplicityAnalysisBase.o: /usr/include/sys/cdefs.h
src/MultiplicityAnalysisBase.o: /usr/include/bits/wordsize.h
src/MultiplicityAnalysisBase.o: /usr/include/gnu/stubs.h
src/MultiplicityAnalysisBase.o: /usr/include/gnu/stubs-64.h
src/MultiplicityAnalysisBase.o: /usr/include/bits/types.h
src/MultiplicityAnalysisBase.o: /usr/include/bits/typesizes.h
src/MultiplicityAnalysisBase.o: /usr/include/libio.h /usr/include/_G_config.h
src/MultiplicityAnalysisBase.o: /usr/include/wchar.h
src/MultiplicityAnalysisBase.o: /usr/include/bits/wchar.h
src/MultiplicityAnalysisBase.o: /usr/include/gconv.h
src/MultiplicityAnalysisBase.o: /usr/include/bits/stdio_lim.h
src/MultiplicityAnalysisBase.o: /usr/include/bits/sys_errlist.h
src/MultiplicityAnalysisBase.o: /usr/include/stdlib.h
src/MultiplicityAnalysisBase.o: /usr/include/sys/types.h /usr/include/time.h
src/MultiplicityAnalysisBase.o: /usr/include/endian.h
src/MultiplicityAnalysisBase.o: /usr/include/bits/endian.h
src/MultiplicityAnalysisBase.o: /usr/include/sys/select.h
src/MultiplicityAnalysisBase.o: /usr/include/bits/select.h
src/MultiplicityAnalysisBase.o: /usr/include/bits/sigset.h
src/MultiplicityAnalysisBase.o: /usr/include/bits/time.h
src/MultiplicityAnalysisBase.o: /usr/include/sys/sysmacros.h
src/MultiplicityAnalysisBase.o: /usr/include/bits/pthreadtypes.h
src/MultiplicityAnalysisBase.o: /usr/include/alloca.h
src/MultiplicityAnalysisBase.o: ./include/MultiplicityAnalysisBase.hh
src/MultiplicityAnalysisBase.o: ./include/base/TreeReader.hh
src/MultiplicityAnalysisBase.o: ./include/base/TreeClassBase.h
src/MultiplicityAnalysisBase.o: ./include/base/UserAnalysisBase.hh
src/MultiplicityAnalysisBase.o: ./include/base/TreeReader.hh
src/MultiplicityAnalysisBase.o: ./include/helper/pdgparticle.hh
src/MultiplicityAnalysis.o: ./include/MultiplicityAnalysis.hh
src/MultiplicityAnalysis.o: ./include/base/TreeReader.hh
src/MultiplicityAnalysis.o: ./include/base/TreeClassBase.h
src/MultiplicityAnalysis.o: ./include/MultiplicityAnalysisBase.hh
src/MultiplicityAnalysis.o: ./include/base/UserAnalysisBase.hh
src/MultiplicityAnalysis.o: ./include/base/TreeReader.hh
src/MultiplicityAnalysis.o: ./include/helper/pdgparticle.hh
src/MultiplicityAnalysis.o: /usr/include/stdlib.h /usr/include/features.h
src/MultiplicityAnalysis.o: /usr/include/sys/cdefs.h
src/MultiplicityAnalysis.o: /usr/include/bits/wordsize.h
src/MultiplicityAnalysis.o: /usr/include/gnu/stubs.h
src/MultiplicityAnalysis.o: /usr/include/gnu/stubs-64.h
src/MultiplicityAnalysis.o: /usr/include/sys/types.h
src/MultiplicityAnalysis.o: /usr/include/bits/types.h
src/MultiplicityAnalysis.o: /usr/include/bits/typesizes.h /usr/include/time.h
src/MultiplicityAnalysis.o: /usr/include/endian.h /usr/include/bits/endian.h
src/MultiplicityAnalysis.o: /usr/include/sys/select.h
src/MultiplicityAnalysis.o: /usr/include/bits/select.h
src/MultiplicityAnalysis.o: /usr/include/bits/sigset.h
src/MultiplicityAnalysis.o: /usr/include/bits/time.h
src/MultiplicityAnalysis.o: /usr/include/sys/sysmacros.h
src/MultiplicityAnalysis.o: /usr/include/bits/pthreadtypes.h
src/MultiplicityAnalysis.o: /usr/include/alloca.h
src/MultiplicityAnalysis.o: ./include/helper/Utilities.hh
src/MultiplicityAnalysis.o: /usr/include/stdio.h /usr/include/libio.h
src/MultiplicityAnalysis.o: /usr/include/_G_config.h /usr/include/wchar.h
src/MultiplicityAnalysis.o: /usr/include/bits/wchar.h /usr/include/gconv.h
src/MultiplicityAnalysis.o: /usr/include/bits/stdio_lim.h
src/MultiplicityAnalysis.o: /usr/include/bits/sys_errlist.h
src/MultiplicityAnalysis.o: ./include/helper/LeptJetStat.h
src/SignificanceAnalysis.o: ./include/SignificanceAnalysis.hh
src/SignificanceAnalysis.o: ./include/base/TreeReader.hh
src/SignificanceAnalysis.o: ./include/base/TreeClassBase.h
src/SignificanceAnalysis.o: ./include/base/UserAnalysisBase.hh
src/SignificanceAnalysis.o: ./include/base/TreeReader.hh
src/SignificanceAnalysis.o: ./include/helper/pdgparticle.hh
src/SignificanceAnalysis.o: /usr/include/stdlib.h /usr/include/features.h
src/SignificanceAnalysis.o: /usr/include/sys/cdefs.h
src/SignificanceAnalysis.o: /usr/include/bits/wordsize.h
src/SignificanceAnalysis.o: /usr/include/gnu/stubs.h
src/SignificanceAnalysis.o: /usr/include/gnu/stubs-64.h
src/SignificanceAnalysis.o: /usr/include/sys/types.h
src/SignificanceAnalysis.o: /usr/include/bits/types.h
src/SignificanceAnalysis.o: /usr/include/bits/typesizes.h /usr/include/time.h
src/SignificanceAnalysis.o: /usr/include/endian.h /usr/include/bits/endian.h
src/SignificanceAnalysis.o: /usr/include/sys/select.h
src/SignificanceAnalysis.o: /usr/include/bits/select.h
src/SignificanceAnalysis.o: /usr/include/bits/sigset.h
src/SignificanceAnalysis.o: /usr/include/bits/time.h
src/SignificanceAnalysis.o: /usr/include/sys/sysmacros.h
src/SignificanceAnalysis.o: /usr/include/bits/pthreadtypes.h
src/SignificanceAnalysis.o: /usr/include/alloca.h
src/SignificanceAnalysis.o: ./include/helper/Utilities.hh
src/SignificanceAnalysis.o: /usr/include/stdio.h /usr/include/libio.h
src/SignificanceAnalysis.o: /usr/include/_G_config.h /usr/include/wchar.h
src/SignificanceAnalysis.o: /usr/include/bits/wchar.h /usr/include/gconv.h
src/SignificanceAnalysis.o: /usr/include/bits/stdio_lim.h
src/SignificanceAnalysis.o: /usr/include/bits/sys_errlist.h
src/PhysQCAnalysis.o: ./include/base/TreeReader.hh
src/PhysQCAnalysis.o: ./include/base/TreeClassBase.h
src/PhysQCAnalysis.o: ./include/helper/AnaClass.hh /usr/include/math.h
src/PhysQCAnalysis.o: /usr/include/features.h /usr/include/sys/cdefs.h
src/PhysQCAnalysis.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
src/PhysQCAnalysis.o: /usr/include/gnu/stubs-64.h
src/PhysQCAnalysis.o: /usr/include/bits/huge_val.h
src/PhysQCAnalysis.o: /usr/include/bits/mathdef.h
src/PhysQCAnalysis.o: /usr/include/bits/mathcalls.h /usr/include/stdlib.h
src/PhysQCAnalysis.o: /usr/include/sys/types.h /usr/include/bits/types.h
src/PhysQCAnalysis.o: /usr/include/bits/typesizes.h /usr/include/time.h
src/PhysQCAnalysis.o: /usr/include/endian.h /usr/include/bits/endian.h
src/PhysQCAnalysis.o: /usr/include/sys/select.h /usr/include/bits/select.h
src/PhysQCAnalysis.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
src/PhysQCAnalysis.o: /usr/include/sys/sysmacros.h
src/PhysQCAnalysis.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
src/PhysQCAnalysis.o: /usr/include/stdio.h /usr/include/libio.h
src/PhysQCAnalysis.o: /usr/include/_G_config.h /usr/include/wchar.h
src/PhysQCAnalysis.o: /usr/include/bits/wchar.h /usr/include/gconv.h
src/PhysQCAnalysis.o: /usr/include/bits/stdio_lim.h
src/PhysQCAnalysis.o: /usr/include/bits/sys_errlist.h
src/PhysQCAnalysis.o: ./include/helper/Utilities.hh
src/PhysQCAnalysis.o: ./include/base/TreeAnalyzerBase.hh
src/PhysQCAnalysis.o: ./include/base/TreeReader.hh
src/PhysQCAnalysis.o: ./include/helper/Utilities.hh ./include/TreeCleaner.hh
src/PhysQCAnalysis.o: ./include/base/UserAnalysisBase.hh
src/PhysQCAnalysis.o: ./include/helper/pdgparticle.hh
src/PhysQCAnalysis.o: ./include/helper/Davismt2.h ./include/PhysQCAnalysis.hh
src/RatioAnalysis.o: ./include/RatioAnalysis.hh ./include/base/TreeReader.hh
src/RatioAnalysis.o: ./include/base/TreeClassBase.h
src/RatioAnalysis.o: ./include/MultiplicityAnalysisBase.hh
src/RatioAnalysis.o: ./include/base/UserAnalysisBase.hh
src/RatioAnalysis.o: ./include/base/TreeReader.hh
src/RatioAnalysis.o: ./include/helper/pdgparticle.hh /usr/include/stdlib.h
src/RatioAnalysis.o: /usr/include/features.h /usr/include/sys/cdefs.h
src/RatioAnalysis.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
src/RatioAnalysis.o: /usr/include/gnu/stubs-64.h /usr/include/sys/types.h
src/RatioAnalysis.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
src/RatioAnalysis.o: /usr/include/time.h /usr/include/endian.h
src/RatioAnalysis.o: /usr/include/bits/endian.h /usr/include/sys/select.h
src/RatioAnalysis.o: /usr/include/bits/select.h /usr/include/bits/sigset.h
src/RatioAnalysis.o: /usr/include/bits/time.h /usr/include/sys/sysmacros.h
src/RatioAnalysis.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
src/RatioAnalysis.o: ./include/helper/Utilities.hh /usr/include/stdio.h
src/RatioAnalysis.o: /usr/include/libio.h /usr/include/_G_config.h
src/RatioAnalysis.o: /usr/include/wchar.h /usr/include/bits/wchar.h
src/RatioAnalysis.o: /usr/include/gconv.h /usr/include/bits/stdio_lim.h
src/RatioAnalysis.o: /usr/include/bits/sys_errlist.h
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
src/LeptJetMultAnalyzer.o: ./include/LeptJetMultAnalyzer.hh
src/LeptJetMultAnalyzer.o: ./include/base/TreeAnalyzerBase.hh
src/LeptJetMultAnalyzer.o: ./include/base/TreeReader.hh
src/LeptJetMultAnalyzer.o: ./include/base/TreeClassBase.h
src/LeptJetMultAnalyzer.o: ./include/helper/Utilities.hh /usr/include/stdio.h
src/LeptJetMultAnalyzer.o: /usr/include/features.h /usr/include/sys/cdefs.h
src/LeptJetMultAnalyzer.o: /usr/include/bits/wordsize.h
src/LeptJetMultAnalyzer.o: /usr/include/gnu/stubs.h
src/LeptJetMultAnalyzer.o: /usr/include/gnu/stubs-64.h
src/LeptJetMultAnalyzer.o: /usr/include/bits/types.h
src/LeptJetMultAnalyzer.o: /usr/include/bits/typesizes.h /usr/include/libio.h
src/LeptJetMultAnalyzer.o: /usr/include/_G_config.h /usr/include/wchar.h
src/LeptJetMultAnalyzer.o: /usr/include/bits/wchar.h /usr/include/gconv.h
src/LeptJetMultAnalyzer.o: /usr/include/bits/stdio_lim.h
src/LeptJetMultAnalyzer.o: /usr/include/bits/sys_errlist.h
src/LeptJetMultAnalyzer.o: /usr/include/stdlib.h /usr/include/sys/types.h
src/LeptJetMultAnalyzer.o: /usr/include/time.h /usr/include/endian.h
src/LeptJetMultAnalyzer.o: /usr/include/bits/endian.h
src/LeptJetMultAnalyzer.o: /usr/include/sys/select.h
src/LeptJetMultAnalyzer.o: /usr/include/bits/select.h
src/LeptJetMultAnalyzer.o: /usr/include/bits/sigset.h
src/LeptJetMultAnalyzer.o: /usr/include/bits/time.h
src/LeptJetMultAnalyzer.o: /usr/include/sys/sysmacros.h
src/LeptJetMultAnalyzer.o: /usr/include/bits/pthreadtypes.h
src/LeptJetMultAnalyzer.o: /usr/include/alloca.h ./include/base/TreeReader.hh
src/LeptJetMultAnalyzer.o: ./include/MultiplicityAnalysis.hh
src/LeptJetMultAnalyzer.o: ./include/MultiplicityAnalysisBase.hh
src/LeptJetMultAnalyzer.o: ./include/base/UserAnalysisBase.hh
src/LeptJetMultAnalyzer.o: ./include/helper/pdgparticle.hh
src/LeptJetMultAnalyzer.o: ./include/helper/LeptJetStat.h
src/LeptJetMultAnalyzer.o: ./include/MassAnalysis.hh
src/LeptJetMultAnalyzer.o: ./include/helper/Davismt2.h /usr/include/math.h
src/LeptJetMultAnalyzer.o: /usr/include/bits/huge_val.h
src/LeptJetMultAnalyzer.o: /usr/include/bits/mathdef.h
src/LeptJetMultAnalyzer.o: /usr/include/bits/mathcalls.h
src/LeptJetMultAnalyzer.o: ./include/helper/Hemisphere.hh
src/LeptJetMultAnalyzer.o: ./include/RatioAnalysis.hh
src/MuonFakeRatioFiller.o: ./include/MuonFakeRatioFiller.hh
src/MuonFakeRatioFiller.o: ./include/base/TreeReader.hh
src/MuonFakeRatioFiller.o: ./include/base/TreeClassBase.h
src/MuonFakeRatioFiller.o: ./include/base/UserAnalysisBase.hh
src/MuonFakeRatioFiller.o: ./include/base/TreeReader.hh
src/MuonFakeRatioFiller.o: ./include/helper/pdgparticle.hh
src/MuonFakeRatioFiller.o: /usr/include/stdlib.h /usr/include/features.h
src/MuonFakeRatioFiller.o: /usr/include/sys/cdefs.h
src/MuonFakeRatioFiller.o: /usr/include/bits/wordsize.h
src/MuonFakeRatioFiller.o: /usr/include/gnu/stubs.h
src/MuonFakeRatioFiller.o: /usr/include/gnu/stubs-64.h
src/MuonFakeRatioFiller.o: /usr/include/sys/types.h /usr/include/bits/types.h
src/MuonFakeRatioFiller.o: /usr/include/bits/typesizes.h /usr/include/time.h
src/MuonFakeRatioFiller.o: /usr/include/endian.h /usr/include/bits/endian.h
src/MuonFakeRatioFiller.o: /usr/include/sys/select.h
src/MuonFakeRatioFiller.o: /usr/include/bits/select.h
src/MuonFakeRatioFiller.o: /usr/include/bits/sigset.h
src/MuonFakeRatioFiller.o: /usr/include/bits/time.h
src/MuonFakeRatioFiller.o: /usr/include/sys/sysmacros.h
src/MuonFakeRatioFiller.o: /usr/include/bits/pthreadtypes.h
src/MuonFakeRatioFiller.o: /usr/include/alloca.h
src/MuonFakeRatioFiller.o: ./include/helper/Utilities.hh /usr/include/stdio.h
src/MuonFakeRatioFiller.o: /usr/include/libio.h /usr/include/_G_config.h
src/MuonFakeRatioFiller.o: /usr/include/wchar.h /usr/include/bits/wchar.h
src/MuonFakeRatioFiller.o: /usr/include/gconv.h /usr/include/bits/stdio_lim.h
src/MuonFakeRatioFiller.o: /usr/include/bits/sys_errlist.h
src/MuonFakeAnalyzer.o: ./include/MuonFakeAnalyzer.hh
src/MuonFakeAnalyzer.o: ./include/base/TreeAnalyzerBase.hh
src/MuonFakeAnalyzer.o: ./include/base/TreeReader.hh
src/MuonFakeAnalyzer.o: ./include/base/TreeClassBase.h
src/MuonFakeAnalyzer.o: ./include/helper/Utilities.hh /usr/include/stdio.h
src/MuonFakeAnalyzer.o: /usr/include/features.h /usr/include/sys/cdefs.h
src/MuonFakeAnalyzer.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
src/MuonFakeAnalyzer.o: /usr/include/gnu/stubs-64.h /usr/include/bits/types.h
src/MuonFakeAnalyzer.o: /usr/include/bits/typesizes.h /usr/include/libio.h
src/MuonFakeAnalyzer.o: /usr/include/_G_config.h /usr/include/wchar.h
src/MuonFakeAnalyzer.o: /usr/include/bits/wchar.h /usr/include/gconv.h
src/MuonFakeAnalyzer.o: /usr/include/bits/stdio_lim.h
src/MuonFakeAnalyzer.o: /usr/include/bits/sys_errlist.h /usr/include/stdlib.h
src/MuonFakeAnalyzer.o: /usr/include/sys/types.h /usr/include/time.h
src/MuonFakeAnalyzer.o: /usr/include/endian.h /usr/include/bits/endian.h
src/MuonFakeAnalyzer.o: /usr/include/sys/select.h /usr/include/bits/select.h
src/MuonFakeAnalyzer.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
src/MuonFakeAnalyzer.o: /usr/include/sys/sysmacros.h
src/MuonFakeAnalyzer.o: /usr/include/bits/pthreadtypes.h
src/MuonFakeAnalyzer.o: /usr/include/alloca.h ./include/base/TreeReader.hh
src/MuonFakeAnalyzer.o: ./include/TreeCleaner.hh
src/MuonFakeAnalyzer.o: ./include/base/UserAnalysisBase.hh
src/MuonFakeAnalyzer.o: ./include/helper/pdgparticle.hh
src/MuonFakeAnalyzer.o: ./include/helper/Davismt2.h /usr/include/math.h
src/MuonFakeAnalyzer.o: /usr/include/bits/huge_val.h
src/MuonFakeAnalyzer.o: /usr/include/bits/mathdef.h
src/MuonFakeAnalyzer.o: /usr/include/bits/mathcalls.h
src/MuonFakeAnalyzer.o: ./include/MuonFakeRatioFiller.hh
src/MuonPlotter.o: ./include/MuonPlotter.hh ./include/helper/AnaClass.hh
src/MuonPlotter.o: /usr/include/math.h /usr/include/features.h
src/MuonPlotter.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
src/MuonPlotter.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h
src/MuonPlotter.o: /usr/include/bits/huge_val.h /usr/include/bits/mathdef.h
src/MuonPlotter.o: /usr/include/bits/mathcalls.h /usr/include/stdlib.h
src/MuonPlotter.o: /usr/include/sys/types.h /usr/include/bits/types.h
src/MuonPlotter.o: /usr/include/bits/typesizes.h /usr/include/time.h
src/MuonPlotter.o: /usr/include/endian.h /usr/include/bits/endian.h
src/MuonPlotter.o: /usr/include/sys/select.h /usr/include/bits/select.h
src/MuonPlotter.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
src/MuonPlotter.o: /usr/include/sys/sysmacros.h
src/MuonPlotter.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
src/MuonPlotter.o: /usr/include/stdio.h /usr/include/libio.h
src/MuonPlotter.o: /usr/include/_G_config.h /usr/include/wchar.h
src/MuonPlotter.o: /usr/include/bits/wchar.h /usr/include/gconv.h
src/MuonPlotter.o: /usr/include/bits/stdio_lim.h
src/MuonPlotter.o: /usr/include/bits/sys_errlist.h
src/MuonPlotter.o: ./include/helper/Utilities.hh ./include/helper/FPRatios.hh
src/MuonPlotter.o: ./include/helper/Utilities.hh
src/helper/FPRatios.o: ./include/helper/FPRatios.hh
