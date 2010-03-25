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

NGLIBS         = $(ROOTGLIBS) -lMinuit -lMinuit2
GLIBS          = $(filter-out -lNew, $(NGLIBS))

SRCS           = src/base/TreeClassBase.C src/base/TreeReader.cc src/base/TreeAnalyzerBase.cc src/base/UserAnalysisBase.cc \
                 src/UserAnalyzer.cc src/TreeAnalyzer.cc src/PhysQCAnalyzer.cc src/TreeSkimmer.cc \
                 src/UserAnalysis.cc src/DiLeptonAnalysis.cc src/TreeCleaner.cc src/MultiplicityAnalysis.cc src/SignificanceAnalysis.cc src/PhysQCAnalysis.cc \
                 src/helper/AnaClass.cc src/helper/Davismt2.cc src/helper/LeptJetStat.cc src/helper/Hemisphere.cc

OBJS           = $(patsubst %.C,%.o,$(SRCS:.cc=.o))

.SUFFIXES: .cc,.C,.hh,.h
.PHONY : clean purge all depend PhysQC

# Rules ====================================
all: RunUserAnalyzer RunTreeAnalyzer RunPhysQCAnalyzer RunTreeSkimmer

RunUserAnalyzer: src/exe/RunUserAnalyzer.C $(OBJS)
	$(CXX) $(CXXFLAGS) -ldl $(GLIBS) $(LDFLAGS) -o $@ $^

RunTreeAnalyzer: src/exe/RunTreeAnalyzer.C $(OBJS)
	$(CXX) $(CXXFLAGS) -ldl $(GLIBS) $(LDFLAGS) -o $@ $^

RunPhysQCAnalyzer: src/exe/RunPhysQCAnalyzer.C $(OBJS)
	$(CXX) $(CXXFLAGS) -ldl $(GLIBS) $(LDFLAGS) -o $@ $^

RunTreeSkimmer: src/exe/RunTreeSkimmer.C $(OBJS)
	$(CXX) $(CXXFLAGS) -ldl $(GLIBS) $(LDFLAGS) -o $@ $^

clean:
	$(RM) $(OBJS)
	$(RM) RunUserAnalyzer
	$(RM) RunTreeAnalyzer
	$(RM) RunPhysQCAnalyzer
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
src/base/TreeAnalyzerBase.o: /usr/include/gnu/stubs.h
src/base/TreeAnalyzerBase.o: /usr/lib/gcc/x86_64-redhat-linux/3.4.6/include/stddef.h
src/base/TreeAnalyzerBase.o: ./include/base/TreeAnalyzerBase.hh
src/base/TreeAnalyzerBase.o: ./include/base/TreeReader.hh
src/base/TreeAnalyzerBase.o: ./include/base/TreeClassBase.h
src/base/TreeAnalyzerBase.o: ./include/helper/Utilities.hh
src/base/TreeAnalyzerBase.o: /usr/include/stdio.h /usr/include/bits/types.h
src/base/TreeAnalyzerBase.o: /usr/include/bits/wordsize.h
src/base/TreeAnalyzerBase.o: /usr/include/bits/typesizes.h
src/base/TreeAnalyzerBase.o: /usr/include/libio.h /usr/include/_G_config.h
src/base/TreeAnalyzerBase.o: /usr/include/wchar.h /usr/include/bits/wchar.h
src/base/TreeAnalyzerBase.o: /usr/include/gconv.h
src/base/TreeAnalyzerBase.o: /usr/lib/gcc/x86_64-redhat-linux/3.4.6/include/stdarg.h
src/base/TreeAnalyzerBase.o: /usr/include/bits/stdio_lim.h
src/base/TreeAnalyzerBase.o: /usr/include/bits/sys_errlist.h
src/base/UserAnalysisBase.o: ./include/base/TreeReader.hh
src/base/UserAnalysisBase.o: ./include/base/TreeClassBase.h
src/base/UserAnalysisBase.o: /usr/include/stdlib.h /usr/include/features.h
src/base/UserAnalysisBase.o: /usr/include/sys/cdefs.h
src/base/UserAnalysisBase.o: /usr/include/gnu/stubs.h
src/base/UserAnalysisBase.o: /usr/lib/gcc/x86_64-redhat-linux/3.4.6/include/stddef.h
src/base/UserAnalysisBase.o: ./include/base/UserAnalysisBase.hh
src/base/UserAnalysisBase.o: ./include/base/TreeReader.hh
src/base/UserAnalysisBase.o: ./include/helper/Utilities.hh
src/base/UserAnalysisBase.o: /usr/include/stdio.h /usr/include/bits/types.h
src/base/UserAnalysisBase.o: /usr/include/bits/wordsize.h
src/base/UserAnalysisBase.o: /usr/include/bits/typesizes.h
src/base/UserAnalysisBase.o: /usr/include/libio.h /usr/include/_G_config.h
src/base/UserAnalysisBase.o: /usr/include/wchar.h /usr/include/bits/wchar.h
src/base/UserAnalysisBase.o: /usr/include/gconv.h
src/base/UserAnalysisBase.o: /usr/lib/gcc/x86_64-redhat-linux/3.4.6/include/stdarg.h
src/base/UserAnalysisBase.o: /usr/include/bits/stdio_lim.h
src/base/UserAnalysisBase.o: /usr/include/bits/sys_errlist.h
src/UserAnalyzer.o: ./include/UserAnalyzer.hh
src/UserAnalyzer.o: ./include/base/TreeAnalyzerBase.hh
src/UserAnalyzer.o: ./include/base/TreeReader.hh
src/UserAnalyzer.o: ./include/base/TreeClassBase.h
src/UserAnalyzer.o: ./include/helper/Utilities.hh /usr/include/stdio.h
src/UserAnalyzer.o: /usr/include/features.h /usr/include/sys/cdefs.h
src/UserAnalyzer.o: /usr/include/gnu/stubs.h
src/UserAnalyzer.o: /usr/lib/gcc/x86_64-redhat-linux/3.4.6/include/stddef.h
src/UserAnalyzer.o: /usr/include/bits/types.h /usr/include/bits/wordsize.h
src/UserAnalyzer.o: /usr/include/bits/typesizes.h /usr/include/libio.h
src/UserAnalyzer.o: /usr/include/_G_config.h /usr/include/wchar.h
src/UserAnalyzer.o: /usr/include/bits/wchar.h /usr/include/gconv.h
src/UserAnalyzer.o: /usr/lib/gcc/x86_64-redhat-linux/3.4.6/include/stdarg.h
src/UserAnalyzer.o: /usr/include/bits/stdio_lim.h
src/UserAnalyzer.o: /usr/include/bits/sys_errlist.h /usr/include/stdlib.h
src/UserAnalyzer.o: ./include/base/TreeReader.hh ./include/UserAnalysis.hh
src/UserAnalyzer.o: ./include/base/UserAnalysisBase.hh
src/TreeAnalyzer.o: ./include/TreeAnalyzer.hh
src/TreeAnalyzer.o: ./include/base/TreeAnalyzerBase.hh
src/TreeAnalyzer.o: ./include/base/TreeReader.hh
src/TreeAnalyzer.o: ./include/base/TreeClassBase.h
src/TreeAnalyzer.o: ./include/helper/Utilities.hh /usr/include/stdio.h
src/TreeAnalyzer.o: /usr/include/features.h /usr/include/sys/cdefs.h
src/TreeAnalyzer.o: /usr/include/gnu/stubs.h
src/TreeAnalyzer.o: /usr/lib/gcc/x86_64-redhat-linux/3.4.6/include/stddef.h
src/TreeAnalyzer.o: /usr/include/bits/types.h /usr/include/bits/wordsize.h
src/TreeAnalyzer.o: /usr/include/bits/typesizes.h /usr/include/libio.h
src/TreeAnalyzer.o: /usr/include/_G_config.h /usr/include/wchar.h
src/TreeAnalyzer.o: /usr/include/bits/wchar.h /usr/include/gconv.h
src/TreeAnalyzer.o: /usr/lib/gcc/x86_64-redhat-linux/3.4.6/include/stdarg.h
src/TreeAnalyzer.o: /usr/include/bits/stdio_lim.h
src/TreeAnalyzer.o: /usr/include/bits/sys_errlist.h /usr/include/stdlib.h
src/TreeAnalyzer.o: ./include/base/TreeReader.hh ./include/TreeCleaner.hh
src/TreeAnalyzer.o: ./include/base/UserAnalysisBase.hh
src/TreeAnalyzer.o: ./include/helper/Davismt2.h /usr/include/math.h
src/TreeAnalyzer.o: /usr/include/bits/huge_val.h /usr/include/bits/mathdef.h
src/TreeAnalyzer.o: /usr/include/bits/mathcalls.h
src/TreeAnalyzer.o: ./include/DiLeptonAnalysis.hh
src/TreeAnalyzer.o: ./include/MultiplicityAnalysis.hh
src/TreeAnalyzer.o: ./include/helper/LeptJetStat.h
src/TreeAnalyzer.o: ./include/SignificanceAnalysis.hh
src/PhysQCAnalyzer.o: ./include/PhysQCAnalyzer.hh
src/PhysQCAnalyzer.o: ./include/base/TreeAnalyzerBase.hh
src/PhysQCAnalyzer.o: ./include/base/TreeReader.hh
src/PhysQCAnalyzer.o: ./include/base/TreeClassBase.h
src/PhysQCAnalyzer.o: ./include/helper/Utilities.hh /usr/include/stdio.h
src/PhysQCAnalyzer.o: /usr/include/features.h /usr/include/sys/cdefs.h
src/PhysQCAnalyzer.o: /usr/include/gnu/stubs.h
src/PhysQCAnalyzer.o: /usr/lib/gcc/x86_64-redhat-linux/3.4.6/include/stddef.h
src/PhysQCAnalyzer.o: /usr/include/bits/types.h /usr/include/bits/wordsize.h
src/PhysQCAnalyzer.o: /usr/include/bits/typesizes.h /usr/include/libio.h
src/PhysQCAnalyzer.o: /usr/include/_G_config.h /usr/include/wchar.h
src/PhysQCAnalyzer.o: /usr/include/bits/wchar.h /usr/include/gconv.h
src/PhysQCAnalyzer.o: /usr/lib/gcc/x86_64-redhat-linux/3.4.6/include/stdarg.h
src/PhysQCAnalyzer.o: /usr/include/bits/stdio_lim.h
src/PhysQCAnalyzer.o: /usr/include/bits/sys_errlist.h /usr/include/stdlib.h
src/PhysQCAnalyzer.o: ./include/base/TreeReader.hh ./include/TreeCleaner.hh
src/PhysQCAnalyzer.o: ./include/base/UserAnalysisBase.hh
src/PhysQCAnalyzer.o: ./include/helper/Davismt2.h /usr/include/math.h
src/PhysQCAnalyzer.o: /usr/include/bits/huge_val.h
src/PhysQCAnalyzer.o: /usr/include/bits/mathdef.h
src/PhysQCAnalyzer.o: /usr/include/bits/mathcalls.h
src/PhysQCAnalyzer.o: ./include/PhysQCAnalysis.hh /usr/include/time.h
src/PhysQCAnalyzer.o: /usr/include/bits/time.h ./include/helper/AnaClass.hh
src/PhysQCAnalyzer.o: ./include/helper/Utilities.hh
src/PhysQCAnalyzer.o: ./include/MultiplicityAnalysis.hh
src/PhysQCAnalyzer.o: ./include/helper/LeptJetStat.h
src/PhysQCAnalyzer.o: ./include/DiLeptonAnalysis.hh
src/TreeSkimmer.o: ./include/TreeSkimmer.hh
src/TreeSkimmer.o: ./include/base/TreeAnalyzerBase.hh
src/TreeSkimmer.o: ./include/base/TreeReader.hh
src/TreeSkimmer.o: ./include/base/TreeClassBase.h
src/TreeSkimmer.o: ./include/helper/Utilities.hh /usr/include/stdio.h
src/TreeSkimmer.o: /usr/include/features.h /usr/include/sys/cdefs.h
src/TreeSkimmer.o: /usr/include/gnu/stubs.h
src/TreeSkimmer.o: /usr/lib/gcc/x86_64-redhat-linux/3.4.6/include/stddef.h
src/TreeSkimmer.o: /usr/include/bits/types.h /usr/include/bits/wordsize.h
src/TreeSkimmer.o: /usr/include/bits/typesizes.h /usr/include/libio.h
src/TreeSkimmer.o: /usr/include/_G_config.h /usr/include/wchar.h
src/TreeSkimmer.o: /usr/include/bits/wchar.h /usr/include/gconv.h
src/TreeSkimmer.o: /usr/lib/gcc/x86_64-redhat-linux/3.4.6/include/stdarg.h
src/TreeSkimmer.o: /usr/include/bits/stdio_lim.h
src/TreeSkimmer.o: /usr/include/bits/sys_errlist.h /usr/include/stdlib.h
src/TreeSkimmer.o: ./include/base/TreeReader.hh
src/UserAnalysis.o: ./include/helper/Utilities.hh /usr/include/stdio.h
src/UserAnalysis.o: /usr/include/features.h /usr/include/sys/cdefs.h
src/UserAnalysis.o: /usr/include/gnu/stubs.h
src/UserAnalysis.o: /usr/lib/gcc/x86_64-redhat-linux/3.4.6/include/stddef.h
src/UserAnalysis.o: /usr/include/bits/types.h /usr/include/bits/wordsize.h
src/UserAnalysis.o: /usr/include/bits/typesizes.h /usr/include/libio.h
src/UserAnalysis.o: /usr/include/_G_config.h /usr/include/wchar.h
src/UserAnalysis.o: /usr/include/bits/wchar.h /usr/include/gconv.h
src/UserAnalysis.o: /usr/lib/gcc/x86_64-redhat-linux/3.4.6/include/stdarg.h
src/UserAnalysis.o: /usr/include/bits/stdio_lim.h
src/UserAnalysis.o: /usr/include/bits/sys_errlist.h /usr/include/stdlib.h
src/UserAnalysis.o: ./include/UserAnalysis.hh ./include/base/TreeReader.hh
src/UserAnalysis.o: ./include/base/TreeClassBase.h
src/UserAnalysis.o: ./include/base/UserAnalysisBase.hh
src/UserAnalysis.o: ./include/base/TreeReader.hh
src/DiLeptonAnalysis.o: ./include/DiLeptonAnalysis.hh
src/DiLeptonAnalysis.o: ./include/base/TreeReader.hh
src/DiLeptonAnalysis.o: ./include/base/TreeClassBase.h
src/DiLeptonAnalysis.o: ./include/base/UserAnalysisBase.hh
src/DiLeptonAnalysis.o: ./include/base/TreeReader.hh
src/DiLeptonAnalysis.o: ./include/helper/Utilities.hh /usr/include/stdio.h
src/DiLeptonAnalysis.o: /usr/include/features.h /usr/include/sys/cdefs.h
src/DiLeptonAnalysis.o: /usr/include/gnu/stubs.h
src/DiLeptonAnalysis.o: /usr/lib/gcc/x86_64-redhat-linux/3.4.6/include/stddef.h
src/DiLeptonAnalysis.o: /usr/include/bits/types.h
src/DiLeptonAnalysis.o: /usr/include/bits/wordsize.h
src/DiLeptonAnalysis.o: /usr/include/bits/typesizes.h /usr/include/libio.h
src/DiLeptonAnalysis.o: /usr/include/_G_config.h /usr/include/wchar.h
src/DiLeptonAnalysis.o: /usr/include/bits/wchar.h /usr/include/gconv.h
src/DiLeptonAnalysis.o: /usr/lib/gcc/x86_64-redhat-linux/3.4.6/include/stdarg.h
src/DiLeptonAnalysis.o: /usr/include/bits/stdio_lim.h
src/DiLeptonAnalysis.o: /usr/include/bits/sys_errlist.h /usr/include/stdlib.h
src/DiLeptonAnalysis.o: ./include/helper/Davismt2.h /usr/include/math.h
src/DiLeptonAnalysis.o: /usr/include/bits/huge_val.h
src/DiLeptonAnalysis.o: /usr/include/bits/mathdef.h
src/DiLeptonAnalysis.o: /usr/include/bits/mathcalls.h
src/TreeCleaner.o: ./include/TreeCleaner.hh ./include/base/TreeReader.hh
src/TreeCleaner.o: ./include/base/TreeClassBase.h
src/TreeCleaner.o: ./include/base/UserAnalysisBase.hh
src/TreeCleaner.o: ./include/base/TreeReader.hh ./include/helper/Utilities.hh
src/TreeCleaner.o: /usr/include/stdio.h /usr/include/features.h
src/TreeCleaner.o: /usr/include/sys/cdefs.h /usr/include/gnu/stubs.h
src/TreeCleaner.o: /usr/lib/gcc/x86_64-redhat-linux/3.4.6/include/stddef.h
src/TreeCleaner.o: /usr/include/bits/types.h /usr/include/bits/wordsize.h
src/TreeCleaner.o: /usr/include/bits/typesizes.h /usr/include/libio.h
src/TreeCleaner.o: /usr/include/_G_config.h /usr/include/wchar.h
src/TreeCleaner.o: /usr/include/bits/wchar.h /usr/include/gconv.h
src/TreeCleaner.o: /usr/lib/gcc/x86_64-redhat-linux/3.4.6/include/stdarg.h
src/TreeCleaner.o: /usr/include/bits/stdio_lim.h
src/TreeCleaner.o: /usr/include/bits/sys_errlist.h /usr/include/stdlib.h
src/TreeCleaner.o: ./include/helper/Davismt2.h /usr/include/math.h
src/TreeCleaner.o: /usr/include/bits/huge_val.h /usr/include/bits/mathdef.h
src/TreeCleaner.o: /usr/include/bits/mathcalls.h
src/MultiplicityAnalysis.o: ./include/MultiplicityAnalysis.hh
src/MultiplicityAnalysis.o: ./include/base/TreeReader.hh
src/MultiplicityAnalysis.o: ./include/base/TreeClassBase.h
src/MultiplicityAnalysis.o: ./include/base/UserAnalysisBase.hh
src/MultiplicityAnalysis.o: ./include/base/TreeReader.hh
src/MultiplicityAnalysis.o: ./include/helper/Utilities.hh
src/MultiplicityAnalysis.o: /usr/include/stdio.h /usr/include/features.h
src/MultiplicityAnalysis.o: /usr/include/sys/cdefs.h /usr/include/gnu/stubs.h
src/MultiplicityAnalysis.o: /usr/lib/gcc/x86_64-redhat-linux/3.4.6/include/stddef.h
src/MultiplicityAnalysis.o: /usr/include/bits/types.h
src/MultiplicityAnalysis.o: /usr/include/bits/wordsize.h
src/MultiplicityAnalysis.o: /usr/include/bits/typesizes.h
src/MultiplicityAnalysis.o: /usr/include/libio.h /usr/include/_G_config.h
src/MultiplicityAnalysis.o: /usr/include/wchar.h /usr/include/bits/wchar.h
src/MultiplicityAnalysis.o: /usr/include/gconv.h
src/MultiplicityAnalysis.o: /usr/lib/gcc/x86_64-redhat-linux/3.4.6/include/stdarg.h
src/MultiplicityAnalysis.o: /usr/include/bits/stdio_lim.h
src/MultiplicityAnalysis.o: /usr/include/bits/sys_errlist.h
src/MultiplicityAnalysis.o: /usr/include/stdlib.h
src/MultiplicityAnalysis.o: ./include/helper/LeptJetStat.h
src/SignificanceAnalysis.o: ./include/SignificanceAnalysis.hh
src/SignificanceAnalysis.o: ./include/base/TreeReader.hh
src/SignificanceAnalysis.o: ./include/base/TreeClassBase.h
src/SignificanceAnalysis.o: ./include/base/UserAnalysisBase.hh
src/SignificanceAnalysis.o: ./include/base/TreeReader.hh
src/SignificanceAnalysis.o: ./include/helper/Utilities.hh
src/SignificanceAnalysis.o: /usr/include/stdio.h /usr/include/features.h
src/SignificanceAnalysis.o: /usr/include/sys/cdefs.h /usr/include/gnu/stubs.h
src/SignificanceAnalysis.o: /usr/lib/gcc/x86_64-redhat-linux/3.4.6/include/stddef.h
src/SignificanceAnalysis.o: /usr/include/bits/types.h
src/SignificanceAnalysis.o: /usr/include/bits/wordsize.h
src/SignificanceAnalysis.o: /usr/include/bits/typesizes.h
src/SignificanceAnalysis.o: /usr/include/libio.h /usr/include/_G_config.h
src/SignificanceAnalysis.o: /usr/include/wchar.h /usr/include/bits/wchar.h
src/SignificanceAnalysis.o: /usr/include/gconv.h
src/SignificanceAnalysis.o: /usr/lib/gcc/x86_64-redhat-linux/3.4.6/include/stdarg.h
src/SignificanceAnalysis.o: /usr/include/bits/stdio_lim.h
src/SignificanceAnalysis.o: /usr/include/bits/sys_errlist.h
src/SignificanceAnalysis.o: /usr/include/stdlib.h
src/PhysQCAnalysis.o: ./include/PhysQCAnalysis.hh /usr/include/time.h
src/PhysQCAnalysis.o: /usr/include/features.h /usr/include/sys/cdefs.h
src/PhysQCAnalysis.o: /usr/include/gnu/stubs.h
src/PhysQCAnalysis.o: /usr/lib/gcc/x86_64-redhat-linux/3.4.6/include/stddef.h
src/PhysQCAnalysis.o: /usr/include/bits/time.h /usr/include/bits/types.h
src/PhysQCAnalysis.o: /usr/include/bits/wordsize.h
src/PhysQCAnalysis.o: /usr/include/bits/typesizes.h
src/PhysQCAnalysis.o: ./include/base/TreeReader.hh
src/PhysQCAnalysis.o: ./include/base/TreeClassBase.h
src/PhysQCAnalysis.o: ./include/helper/AnaClass.hh /usr/include/math.h
src/PhysQCAnalysis.o: /usr/include/bits/huge_val.h
src/PhysQCAnalysis.o: /usr/include/bits/mathdef.h
src/PhysQCAnalysis.o: /usr/include/bits/mathcalls.h /usr/include/stdlib.h
src/PhysQCAnalysis.o: /usr/include/stdio.h /usr/include/libio.h
src/PhysQCAnalysis.o: /usr/include/_G_config.h /usr/include/wchar.h
src/PhysQCAnalysis.o: /usr/include/bits/wchar.h /usr/include/gconv.h
src/PhysQCAnalysis.o: /usr/lib/gcc/x86_64-redhat-linux/3.4.6/include/stdarg.h
src/PhysQCAnalysis.o: /usr/include/bits/stdio_lim.h
src/PhysQCAnalysis.o: /usr/include/bits/sys_errlist.h
src/PhysQCAnalysis.o: ./include/helper/Utilities.hh
src/PhysQCAnalysis.o: ./include/base/UserAnalysisBase.hh
src/PhysQCAnalysis.o: ./include/base/TreeReader.hh
src/PhysQCAnalysis.o: ./include/helper/Utilities.hh ./include/TreeCleaner.hh
src/PhysQCAnalysis.o: ./include/helper/Davismt2.h
src/PhysQCAnalysis.o: ./include/base/TreeAnalyzerBase.hh
src/helper/AnaClass.o: ./include/helper/AnaClass.hh /usr/include/math.h
src/helper/AnaClass.o: /usr/include/features.h /usr/include/sys/cdefs.h
src/helper/AnaClass.o: /usr/include/gnu/stubs.h /usr/include/bits/huge_val.h
src/helper/AnaClass.o: /usr/include/bits/mathdef.h
src/helper/AnaClass.o: /usr/include/bits/mathcalls.h /usr/include/stdlib.h
src/helper/AnaClass.o: /usr/lib/gcc/x86_64-redhat-linux/3.4.6/include/stddef.h
src/helper/AnaClass.o: /usr/include/stdio.h /usr/include/bits/types.h
src/helper/AnaClass.o: /usr/include/bits/wordsize.h
src/helper/AnaClass.o: /usr/include/bits/typesizes.h /usr/include/libio.h
src/helper/AnaClass.o: /usr/include/_G_config.h /usr/include/wchar.h
src/helper/AnaClass.o: /usr/include/bits/wchar.h /usr/include/gconv.h
src/helper/AnaClass.o: /usr/lib/gcc/x86_64-redhat-linux/3.4.6/include/stdarg.h
src/helper/AnaClass.o: /usr/include/bits/stdio_lim.h
src/helper/AnaClass.o: /usr/include/bits/sys_errlist.h /usr/include/time.h
src/helper/AnaClass.o: /usr/include/bits/time.h ./include/helper/Utilities.hh
src/helper/AnaClass.o: ./include/helper/Utilities.hh
src/helper/Davismt2.o: ./include/helper/Davismt2.h /usr/include/math.h
src/helper/Davismt2.o: /usr/include/features.h /usr/include/sys/cdefs.h
src/helper/Davismt2.o: /usr/include/gnu/stubs.h /usr/include/bits/huge_val.h
src/helper/Davismt2.o: /usr/include/bits/mathdef.h
src/helper/Davismt2.o: /usr/include/bits/mathcalls.h
src/helper/LeptJetStat.o: ./include/helper/LeptJetStat.h
src/helper/Hemisphere.o: ./include/helper/Hemisphere.hh
src/helper/Hemisphere.o: ./include/helper/Utilities.hh /usr/include/stdio.h
src/helper/Hemisphere.o: /usr/include/features.h /usr/include/sys/cdefs.h
src/helper/Hemisphere.o: /usr/include/gnu/stubs.h
src/helper/Hemisphere.o: /usr/lib/gcc/x86_64-redhat-linux/3.4.6/include/stddef.h
src/helper/Hemisphere.o: /usr/include/bits/types.h
src/helper/Hemisphere.o: /usr/include/bits/wordsize.h
src/helper/Hemisphere.o: /usr/include/bits/typesizes.h /usr/include/libio.h
src/helper/Hemisphere.o: /usr/include/_G_config.h /usr/include/wchar.h
src/helper/Hemisphere.o: /usr/include/bits/wchar.h /usr/include/gconv.h
src/helper/Hemisphere.o: /usr/lib/gcc/x86_64-redhat-linux/3.4.6/include/stdarg.h
src/helper/Hemisphere.o: /usr/include/bits/stdio_lim.h
src/helper/Hemisphere.o: /usr/include/bits/sys_errlist.h
src/helper/Hemisphere.o: /usr/include/stdlib.h
