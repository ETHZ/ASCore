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
                 src/helper/TMctLib.cc src/helper/mctlib.cc src/helper/FPRatios.cc \
                 src/helper/AnaClass.cc src/helper/Davismt2.cc src/helper/LeptJetStat.cc src/helper/Hemisphere.cc  src/LeptJetMultAnalyzer.cc \
                 src/MuonPlotter.cc src/MassPlotter.cc src/helper/MetaTreeClassBase.C src/helper/MassAnalysisTreeClass.C \
                 src/JZBAnalyzer.cc src/JZBAnalysis.cc \
                 src/SSDLAnalyzer.cc src/SSDLAnalysis.cc

OBJS           = $(patsubst %.C,%.o,$(SRCS:.cc=.o))

.SUFFIXES: .cc,.C,.hh,.h
.PHONY : clean purge all depend PhysQC

# Rules ====================================
all: RunUserAnalyzer RunTreeAnalyzer RunPhysQCAnalyzer RunTreeSkimmer RunLeptJetMultAnalyzer MakeMuonPlots RunJZBAnalyzer RunSSDLAnalyzer MakeMassPlots

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

src/base/TreeClassBase.o: ./include/base/TreeClassBase.h
src/base/TreeReader.o: ./include/base/TreeReader.hh
src/base/TreeReader.o: ./include/base/TreeClassBase.h
src/base/TreeAnalyzerBase.o: /usr/include/stdlib.h
src/base/TreeAnalyzerBase.o: /usr/include/Availability.h
src/base/TreeAnalyzerBase.o: /usr/include/AvailabilityInternal.h
src/base/TreeAnalyzerBase.o: /usr/include/_types.h /usr/include/sys/_types.h
src/base/TreeAnalyzerBase.o: /usr/include/sys/cdefs.h
src/base/TreeAnalyzerBase.o: /usr/include/machine/_types.h
src/base/TreeAnalyzerBase.o: /usr/include/i386/_types.h
src/base/TreeAnalyzerBase.o: /usr/include/sys/wait.h
src/base/TreeAnalyzerBase.o: /usr/include/sys/signal.h
src/base/TreeAnalyzerBase.o: /usr/include/sys/appleapiopts.h
src/base/TreeAnalyzerBase.o: /usr/include/machine/signal.h
src/base/TreeAnalyzerBase.o: /usr/include/i386/signal.h
src/base/TreeAnalyzerBase.o: /usr/include/i386/_structs.h
src/base/TreeAnalyzerBase.o: /usr/include/sys/_structs.h
src/base/TreeAnalyzerBase.o: /usr/include/machine/_structs.h
src/base/TreeAnalyzerBase.o: /usr/include/sys/resource.h
src/base/TreeAnalyzerBase.o: /usr/include/machine/endian.h
src/base/TreeAnalyzerBase.o: /usr/include/i386/endian.h
src/base/TreeAnalyzerBase.o: /usr/include/sys/_endian.h
src/base/TreeAnalyzerBase.o: /usr/include/libkern/_OSByteOrder.h
src/base/TreeAnalyzerBase.o: /usr/include/libkern/i386/_OSByteOrder.h
src/base/TreeAnalyzerBase.o: /usr/include/alloca.h
src/base/TreeAnalyzerBase.o: /usr/include/machine/types.h
src/base/TreeAnalyzerBase.o: /usr/include/i386/types.h
src/base/TreeAnalyzerBase.o: ./include/base/TreeAnalyzerBase.hh
src/base/TreeAnalyzerBase.o: ./include/base/TreeReader.hh
src/base/TreeAnalyzerBase.o: ./include/base/TreeClassBase.h
src/base/TreeAnalyzerBase.o: ./include/helper/Utilities.hh
src/base/TreeAnalyzerBase.o: /usr/include/stdio.h
src/base/TreeAnalyzerBase.o: /usr/include/secure/_stdio.h
src/base/TreeAnalyzerBase.o: /usr/include/secure/_common.h
src/base/UserAnalysisBase.o: ./include/base/TreeReader.hh
src/base/UserAnalysisBase.o: ./include/base/TreeClassBase.h
src/base/UserAnalysisBase.o: /usr/include/stdlib.h
src/base/UserAnalysisBase.o: /usr/include/Availability.h
src/base/UserAnalysisBase.o: /usr/include/AvailabilityInternal.h
src/base/UserAnalysisBase.o: /usr/include/_types.h /usr/include/sys/_types.h
src/base/UserAnalysisBase.o: /usr/include/sys/cdefs.h
src/base/UserAnalysisBase.o: /usr/include/machine/_types.h
src/base/UserAnalysisBase.o: /usr/include/i386/_types.h
src/base/UserAnalysisBase.o: /usr/include/sys/wait.h
src/base/UserAnalysisBase.o: /usr/include/sys/signal.h
src/base/UserAnalysisBase.o: /usr/include/sys/appleapiopts.h
src/base/UserAnalysisBase.o: /usr/include/machine/signal.h
src/base/UserAnalysisBase.o: /usr/include/i386/signal.h
src/base/UserAnalysisBase.o: /usr/include/i386/_structs.h
src/base/UserAnalysisBase.o: /usr/include/sys/_structs.h
src/base/UserAnalysisBase.o: /usr/include/machine/_structs.h
src/base/UserAnalysisBase.o: /usr/include/sys/resource.h
src/base/UserAnalysisBase.o: /usr/include/machine/endian.h
src/base/UserAnalysisBase.o: /usr/include/i386/endian.h
src/base/UserAnalysisBase.o: /usr/include/sys/_endian.h
src/base/UserAnalysisBase.o: /usr/include/libkern/_OSByteOrder.h
src/base/UserAnalysisBase.o: /usr/include/libkern/i386/_OSByteOrder.h
src/base/UserAnalysisBase.o: /usr/include/alloca.h
src/base/UserAnalysisBase.o: /usr/include/machine/types.h
src/base/UserAnalysisBase.o: /usr/include/i386/types.h
src/base/UserAnalysisBase.o: ./include/helper/pdgparticle.hh
src/base/UserAnalysisBase.o: ./include/base/UserAnalysisBase.hh
src/base/UserAnalysisBase.o: ./include/base/TreeReader.hh
src/base/UserAnalysisBase.o: ./include/helper/Utilities.hh
src/base/UserAnalysisBase.o: /usr/include/stdio.h
src/base/UserAnalysisBase.o: /usr/include/secure/_stdio.h
src/base/UserAnalysisBase.o: /usr/include/secure/_common.h
src/UserAnalyzer.o: ./include/UserAnalyzer.hh
src/UserAnalyzer.o: ./include/base/TreeAnalyzerBase.hh
src/UserAnalyzer.o: ./include/base/TreeReader.hh
src/UserAnalyzer.o: ./include/base/TreeClassBase.h
src/UserAnalyzer.o: ./include/helper/Utilities.hh /usr/include/stdio.h
src/UserAnalyzer.o: /usr/include/_types.h /usr/include/sys/_types.h
src/UserAnalyzer.o: /usr/include/sys/cdefs.h /usr/include/machine/_types.h
src/UserAnalyzer.o: /usr/include/i386/_types.h /usr/include/secure/_stdio.h
src/UserAnalyzer.o: /usr/include/secure/_common.h /usr/include/stdlib.h
src/UserAnalyzer.o: /usr/include/Availability.h
src/UserAnalyzer.o: /usr/include/AvailabilityInternal.h
src/UserAnalyzer.o: /usr/include/sys/wait.h /usr/include/sys/signal.h
src/UserAnalyzer.o: /usr/include/sys/appleapiopts.h
src/UserAnalyzer.o: /usr/include/machine/signal.h /usr/include/i386/signal.h
src/UserAnalyzer.o: /usr/include/i386/_structs.h /usr/include/sys/_structs.h
src/UserAnalyzer.o: /usr/include/machine/_structs.h
src/UserAnalyzer.o: /usr/include/sys/resource.h /usr/include/machine/endian.h
src/UserAnalyzer.o: /usr/include/i386/endian.h /usr/include/sys/_endian.h
src/UserAnalyzer.o: /usr/include/libkern/_OSByteOrder.h
src/UserAnalyzer.o: /usr/include/libkern/i386/_OSByteOrder.h
src/UserAnalyzer.o: /usr/include/alloca.h /usr/include/machine/types.h
src/UserAnalyzer.o: /usr/include/i386/types.h ./include/base/TreeReader.hh
src/UserAnalyzer.o: ./include/UserAnalysis.hh
src/UserAnalyzer.o: ./include/base/UserAnalysisBase.hh
src/UserAnalyzer.o: ./include/helper/pdgparticle.hh
src/TreeAnalyzer.o: ./include/TreeAnalyzer.hh
src/TreeAnalyzer.o: ./include/base/TreeAnalyzerBase.hh
src/TreeAnalyzer.o: ./include/base/TreeReader.hh
src/TreeAnalyzer.o: ./include/base/TreeClassBase.h
src/TreeAnalyzer.o: ./include/helper/Utilities.hh /usr/include/stdio.h
src/TreeAnalyzer.o: /usr/include/_types.h /usr/include/sys/_types.h
src/TreeAnalyzer.o: /usr/include/sys/cdefs.h /usr/include/machine/_types.h
src/TreeAnalyzer.o: /usr/include/i386/_types.h /usr/include/secure/_stdio.h
src/TreeAnalyzer.o: /usr/include/secure/_common.h /usr/include/stdlib.h
src/TreeAnalyzer.o: /usr/include/Availability.h
src/TreeAnalyzer.o: /usr/include/AvailabilityInternal.h
src/TreeAnalyzer.o: /usr/include/sys/wait.h /usr/include/sys/signal.h
src/TreeAnalyzer.o: /usr/include/sys/appleapiopts.h
src/TreeAnalyzer.o: /usr/include/machine/signal.h /usr/include/i386/signal.h
src/TreeAnalyzer.o: /usr/include/i386/_structs.h /usr/include/sys/_structs.h
src/TreeAnalyzer.o: /usr/include/machine/_structs.h
src/TreeAnalyzer.o: /usr/include/sys/resource.h /usr/include/machine/endian.h
src/TreeAnalyzer.o: /usr/include/i386/endian.h /usr/include/sys/_endian.h
src/TreeAnalyzer.o: /usr/include/libkern/_OSByteOrder.h
src/TreeAnalyzer.o: /usr/include/libkern/i386/_OSByteOrder.h
src/TreeAnalyzer.o: /usr/include/alloca.h /usr/include/machine/types.h
src/TreeAnalyzer.o: /usr/include/i386/types.h ./include/base/TreeReader.hh
src/TreeAnalyzer.o: ./include/TreeCleaner.hh
src/TreeAnalyzer.o: ./include/base/UserAnalysisBase.hh
src/TreeAnalyzer.o: ./include/helper/pdgparticle.hh
src/TreeAnalyzer.o: ./include/helper/Davismt2.h /usr/include/math.h
src/TreeAnalyzer.o: /usr/include/architecture/i386/math.h
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
src/PhysQCAnalyzer.o: /usr/include/_types.h /usr/include/sys/_types.h
src/PhysQCAnalyzer.o: /usr/include/sys/cdefs.h /usr/include/machine/_types.h
src/PhysQCAnalyzer.o: /usr/include/i386/_types.h /usr/include/secure/_stdio.h
src/PhysQCAnalyzer.o: /usr/include/secure/_common.h /usr/include/stdlib.h
src/PhysQCAnalyzer.o: /usr/include/Availability.h
src/PhysQCAnalyzer.o: /usr/include/AvailabilityInternal.h
src/PhysQCAnalyzer.o: /usr/include/sys/wait.h /usr/include/sys/signal.h
src/PhysQCAnalyzer.o: /usr/include/sys/appleapiopts.h
src/PhysQCAnalyzer.o: /usr/include/machine/signal.h
src/PhysQCAnalyzer.o: /usr/include/i386/signal.h /usr/include/i386/_structs.h
src/PhysQCAnalyzer.o: /usr/include/sys/_structs.h
src/PhysQCAnalyzer.o: /usr/include/machine/_structs.h
src/PhysQCAnalyzer.o: /usr/include/sys/resource.h
src/PhysQCAnalyzer.o: /usr/include/machine/endian.h
src/PhysQCAnalyzer.o: /usr/include/i386/endian.h /usr/include/sys/_endian.h
src/PhysQCAnalyzer.o: /usr/include/libkern/_OSByteOrder.h
src/PhysQCAnalyzer.o: /usr/include/libkern/i386/_OSByteOrder.h
src/PhysQCAnalyzer.o: /usr/include/alloca.h /usr/include/machine/types.h
src/PhysQCAnalyzer.o: /usr/include/i386/types.h ./include/base/TreeReader.hh
src/PhysQCAnalyzer.o: ./include/TreeCleaner.hh
src/PhysQCAnalyzer.o: ./include/base/UserAnalysisBase.hh
src/PhysQCAnalyzer.o: ./include/helper/pdgparticle.hh
src/PhysQCAnalyzer.o: ./include/helper/Davismt2.h /usr/include/math.h
src/PhysQCAnalyzer.o: /usr/include/architecture/i386/math.h
src/PhysQCAnalyzer.o: ./include/PhysQCAnalysis.hh /usr/include/time.h
src/PhysQCAnalyzer.o: /usr/include/_structs.h ./include/helper/AnaClass.hh
src/PhysQCAnalyzer.o: ./include/helper/Utilities.hh
src/PhysQCAnalyzer.o: ./include/helper/MetaTreeClassBase.h
src/PhysQCAnalyzer.o: ./include/MultiplicityAnalysis.hh
src/PhysQCAnalyzer.o: ./include/MultiplicityAnalysisBase.hh
src/PhysQCAnalyzer.o: ./include/helper/LeptJetStat.h
src/PhysQCAnalyzer.o: ./include/DiLeptonAnalysis.hh
src/TreeSkimmer.o: ./include/TreeSkimmer.hh
src/TreeSkimmer.o: ./include/base/TreeAnalyzerBase.hh
src/TreeSkimmer.o: ./include/base/TreeReader.hh
src/TreeSkimmer.o: ./include/base/TreeClassBase.h
src/TreeSkimmer.o: ./include/helper/Utilities.hh /usr/include/stdio.h
src/TreeSkimmer.o: /usr/include/_types.h /usr/include/sys/_types.h
src/TreeSkimmer.o: /usr/include/sys/cdefs.h /usr/include/machine/_types.h
src/TreeSkimmer.o: /usr/include/i386/_types.h /usr/include/secure/_stdio.h
src/TreeSkimmer.o: /usr/include/secure/_common.h /usr/include/stdlib.h
src/TreeSkimmer.o: /usr/include/Availability.h
src/TreeSkimmer.o: /usr/include/AvailabilityInternal.h
src/TreeSkimmer.o: /usr/include/sys/wait.h /usr/include/sys/signal.h
src/TreeSkimmer.o: /usr/include/sys/appleapiopts.h
src/TreeSkimmer.o: /usr/include/machine/signal.h /usr/include/i386/signal.h
src/TreeSkimmer.o: /usr/include/i386/_structs.h /usr/include/sys/_structs.h
src/TreeSkimmer.o: /usr/include/machine/_structs.h
src/TreeSkimmer.o: /usr/include/sys/resource.h /usr/include/machine/endian.h
src/TreeSkimmer.o: /usr/include/i386/endian.h /usr/include/sys/_endian.h
src/TreeSkimmer.o: /usr/include/libkern/_OSByteOrder.h
src/TreeSkimmer.o: /usr/include/libkern/i386/_OSByteOrder.h
src/TreeSkimmer.o: /usr/include/alloca.h /usr/include/machine/types.h
src/TreeSkimmer.o: /usr/include/i386/types.h ./include/base/TreeReader.hh
src/MassAnalysis.o: ./include/helper/Utilities.hh /usr/include/stdio.h
src/MassAnalysis.o: /usr/include/_types.h /usr/include/sys/_types.h
src/MassAnalysis.o: /usr/include/sys/cdefs.h /usr/include/machine/_types.h
src/MassAnalysis.o: /usr/include/i386/_types.h /usr/include/secure/_stdio.h
src/MassAnalysis.o: /usr/include/secure/_common.h /usr/include/stdlib.h
src/MassAnalysis.o: /usr/include/Availability.h
src/MassAnalysis.o: /usr/include/AvailabilityInternal.h
src/MassAnalysis.o: /usr/include/sys/wait.h /usr/include/sys/signal.h
src/MassAnalysis.o: /usr/include/sys/appleapiopts.h
src/MassAnalysis.o: /usr/include/machine/signal.h /usr/include/i386/signal.h
src/MassAnalysis.o: /usr/include/i386/_structs.h /usr/include/sys/_structs.h
src/MassAnalysis.o: /usr/include/machine/_structs.h
src/MassAnalysis.o: /usr/include/sys/resource.h /usr/include/machine/endian.h
src/MassAnalysis.o: /usr/include/i386/endian.h /usr/include/sys/_endian.h
src/MassAnalysis.o: /usr/include/libkern/_OSByteOrder.h
src/MassAnalysis.o: /usr/include/libkern/i386/_OSByteOrder.h
src/MassAnalysis.o: /usr/include/alloca.h /usr/include/machine/types.h
src/MassAnalysis.o: /usr/include/i386/types.h ./include/MassAnalysis.hh
src/MassAnalysis.o: ./include/base/TreeReader.hh
src/MassAnalysis.o: ./include/base/TreeClassBase.h
src/MassAnalysis.o: ./include/MultiplicityAnalysisBase.hh
src/MassAnalysis.o: ./include/base/UserAnalysisBase.hh
src/MassAnalysis.o: ./include/base/TreeReader.hh
src/MassAnalysis.o: ./include/helper/pdgparticle.hh
src/MassAnalysis.o: ./include/helper/Davismt2.h /usr/include/math.h
src/MassAnalysis.o: /usr/include/architecture/i386/math.h
src/MassAnalysis.o: ./include/helper/TMctLib.h ./include/helper/mctlib.h
src/MassAnalysis.o: ./include/helper/Hemisphere.hh
src/UserAnalysis.o: ./include/helper/Utilities.hh /usr/include/stdio.h
src/UserAnalysis.o: /usr/include/_types.h /usr/include/sys/_types.h
src/UserAnalysis.o: /usr/include/sys/cdefs.h /usr/include/machine/_types.h
src/UserAnalysis.o: /usr/include/i386/_types.h /usr/include/secure/_stdio.h
src/UserAnalysis.o: /usr/include/secure/_common.h /usr/include/stdlib.h
src/UserAnalysis.o: /usr/include/Availability.h
src/UserAnalysis.o: /usr/include/AvailabilityInternal.h
src/UserAnalysis.o: /usr/include/sys/wait.h /usr/include/sys/signal.h
src/UserAnalysis.o: /usr/include/sys/appleapiopts.h
src/UserAnalysis.o: /usr/include/machine/signal.h /usr/include/i386/signal.h
src/UserAnalysis.o: /usr/include/i386/_structs.h /usr/include/sys/_structs.h
src/UserAnalysis.o: /usr/include/machine/_structs.h
src/UserAnalysis.o: /usr/include/sys/resource.h /usr/include/machine/endian.h
src/UserAnalysis.o: /usr/include/i386/endian.h /usr/include/sys/_endian.h
src/UserAnalysis.o: /usr/include/libkern/_OSByteOrder.h
src/UserAnalysis.o: /usr/include/libkern/i386/_OSByteOrder.h
src/UserAnalysis.o: /usr/include/alloca.h /usr/include/machine/types.h
src/UserAnalysis.o: /usr/include/i386/types.h ./include/UserAnalysis.hh
src/UserAnalysis.o: ./include/base/TreeReader.hh
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
src/DiLeptonAnalysis.o: /usr/include/Availability.h
src/DiLeptonAnalysis.o: /usr/include/AvailabilityInternal.h
src/DiLeptonAnalysis.o: /usr/include/_types.h /usr/include/sys/_types.h
src/DiLeptonAnalysis.o: /usr/include/sys/cdefs.h
src/DiLeptonAnalysis.o: /usr/include/machine/_types.h
src/DiLeptonAnalysis.o: /usr/include/i386/_types.h /usr/include/sys/wait.h
src/DiLeptonAnalysis.o: /usr/include/sys/signal.h
src/DiLeptonAnalysis.o: /usr/include/sys/appleapiopts.h
src/DiLeptonAnalysis.o: /usr/include/machine/signal.h
src/DiLeptonAnalysis.o: /usr/include/i386/signal.h
src/DiLeptonAnalysis.o: /usr/include/i386/_structs.h
src/DiLeptonAnalysis.o: /usr/include/sys/_structs.h
src/DiLeptonAnalysis.o: /usr/include/machine/_structs.h
src/DiLeptonAnalysis.o: /usr/include/sys/resource.h
src/DiLeptonAnalysis.o: /usr/include/machine/endian.h
src/DiLeptonAnalysis.o: /usr/include/i386/endian.h /usr/include/sys/_endian.h
src/DiLeptonAnalysis.o: /usr/include/libkern/_OSByteOrder.h
src/DiLeptonAnalysis.o: /usr/include/libkern/i386/_OSByteOrder.h
src/DiLeptonAnalysis.o: /usr/include/alloca.h /usr/include/machine/types.h
src/DiLeptonAnalysis.o: /usr/include/i386/types.h
src/DiLeptonAnalysis.o: ./include/helper/Utilities.hh /usr/include/stdio.h
src/DiLeptonAnalysis.o: /usr/include/secure/_stdio.h
src/DiLeptonAnalysis.o: /usr/include/secure/_common.h
src/DiLeptonAnalysis.o: ./include/helper/Davismt2.h /usr/include/math.h
src/DiLeptonAnalysis.o: /usr/include/architecture/i386/math.h
src/TreeCleaner.o: ./include/TreeCleaner.hh ./include/base/TreeReader.hh
src/TreeCleaner.o: ./include/base/TreeClassBase.h
src/TreeCleaner.o: ./include/base/UserAnalysisBase.hh
src/TreeCleaner.o: ./include/base/TreeReader.hh
src/TreeCleaner.o: ./include/helper/pdgparticle.hh /usr/include/stdlib.h
src/TreeCleaner.o: /usr/include/Availability.h
src/TreeCleaner.o: /usr/include/AvailabilityInternal.h /usr/include/_types.h
src/TreeCleaner.o: /usr/include/sys/_types.h /usr/include/sys/cdefs.h
src/TreeCleaner.o: /usr/include/machine/_types.h /usr/include/i386/_types.h
src/TreeCleaner.o: /usr/include/sys/wait.h /usr/include/sys/signal.h
src/TreeCleaner.o: /usr/include/sys/appleapiopts.h
src/TreeCleaner.o: /usr/include/machine/signal.h /usr/include/i386/signal.h
src/TreeCleaner.o: /usr/include/i386/_structs.h /usr/include/sys/_structs.h
src/TreeCleaner.o: /usr/include/machine/_structs.h
src/TreeCleaner.o: /usr/include/sys/resource.h /usr/include/machine/endian.h
src/TreeCleaner.o: /usr/include/i386/endian.h /usr/include/sys/_endian.h
src/TreeCleaner.o: /usr/include/libkern/_OSByteOrder.h
src/TreeCleaner.o: /usr/include/libkern/i386/_OSByteOrder.h
src/TreeCleaner.o: /usr/include/alloca.h /usr/include/machine/types.h
src/TreeCleaner.o: /usr/include/i386/types.h ./include/helper/Utilities.hh
src/TreeCleaner.o: /usr/include/stdio.h /usr/include/secure/_stdio.h
src/TreeCleaner.o: /usr/include/secure/_common.h ./include/helper/Davismt2.h
src/TreeCleaner.o: /usr/include/math.h /usr/include/architecture/i386/math.h
src/MultiplicityAnalysisBase.o: ./include/helper/Utilities.hh
src/MultiplicityAnalysisBase.o: /usr/include/stdio.h /usr/include/_types.h
src/MultiplicityAnalysisBase.o: /usr/include/sys/_types.h
src/MultiplicityAnalysisBase.o: /usr/include/sys/cdefs.h
src/MultiplicityAnalysisBase.o: /usr/include/machine/_types.h
src/MultiplicityAnalysisBase.o: /usr/include/i386/_types.h
src/MultiplicityAnalysisBase.o: /usr/include/secure/_stdio.h
src/MultiplicityAnalysisBase.o: /usr/include/secure/_common.h
src/MultiplicityAnalysisBase.o: /usr/include/stdlib.h
src/MultiplicityAnalysisBase.o: /usr/include/Availability.h
src/MultiplicityAnalysisBase.o: /usr/include/AvailabilityInternal.h
src/MultiplicityAnalysisBase.o: /usr/include/sys/wait.h
src/MultiplicityAnalysisBase.o: /usr/include/sys/signal.h
src/MultiplicityAnalysisBase.o: /usr/include/sys/appleapiopts.h
src/MultiplicityAnalysisBase.o: /usr/include/machine/signal.h
src/MultiplicityAnalysisBase.o: /usr/include/i386/signal.h
src/MultiplicityAnalysisBase.o: /usr/include/i386/_structs.h
src/MultiplicityAnalysisBase.o: /usr/include/sys/_structs.h
src/MultiplicityAnalysisBase.o: /usr/include/machine/_structs.h
src/MultiplicityAnalysisBase.o: /usr/include/sys/resource.h
src/MultiplicityAnalysisBase.o: /usr/include/machine/endian.h
src/MultiplicityAnalysisBase.o: /usr/include/i386/endian.h
src/MultiplicityAnalysisBase.o: /usr/include/sys/_endian.h
src/MultiplicityAnalysisBase.o: /usr/include/libkern/_OSByteOrder.h
src/MultiplicityAnalysisBase.o: /usr/include/libkern/i386/_OSByteOrder.h
src/MultiplicityAnalysisBase.o: /usr/include/alloca.h
src/MultiplicityAnalysisBase.o: /usr/include/machine/types.h
src/MultiplicityAnalysisBase.o: /usr/include/i386/types.h
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
src/MultiplicityAnalysis.o: /usr/include/stdlib.h /usr/include/Availability.h
src/MultiplicityAnalysis.o: /usr/include/AvailabilityInternal.h
src/MultiplicityAnalysis.o: /usr/include/_types.h /usr/include/sys/_types.h
src/MultiplicityAnalysis.o: /usr/include/sys/cdefs.h
src/MultiplicityAnalysis.o: /usr/include/machine/_types.h
src/MultiplicityAnalysis.o: /usr/include/i386/_types.h
src/MultiplicityAnalysis.o: /usr/include/sys/wait.h /usr/include/sys/signal.h
src/MultiplicityAnalysis.o: /usr/include/sys/appleapiopts.h
src/MultiplicityAnalysis.o: /usr/include/machine/signal.h
src/MultiplicityAnalysis.o: /usr/include/i386/signal.h
src/MultiplicityAnalysis.o: /usr/include/i386/_structs.h
src/MultiplicityAnalysis.o: /usr/include/sys/_structs.h
src/MultiplicityAnalysis.o: /usr/include/machine/_structs.h
src/MultiplicityAnalysis.o: /usr/include/sys/resource.h
src/MultiplicityAnalysis.o: /usr/include/machine/endian.h
src/MultiplicityAnalysis.o: /usr/include/i386/endian.h
src/MultiplicityAnalysis.o: /usr/include/sys/_endian.h
src/MultiplicityAnalysis.o: /usr/include/libkern/_OSByteOrder.h
src/MultiplicityAnalysis.o: /usr/include/libkern/i386/_OSByteOrder.h
src/MultiplicityAnalysis.o: /usr/include/alloca.h
src/MultiplicityAnalysis.o: /usr/include/machine/types.h
src/MultiplicityAnalysis.o: /usr/include/i386/types.h
src/MultiplicityAnalysis.o: ./include/helper/Utilities.hh
src/MultiplicityAnalysis.o: /usr/include/stdio.h /usr/include/secure/_stdio.h
src/MultiplicityAnalysis.o: /usr/include/secure/_common.h
src/MultiplicityAnalysis.o: ./include/helper/LeptJetStat.h
src/SignificanceAnalysis.o: ./include/SignificanceAnalysis.hh
src/SignificanceAnalysis.o: ./include/base/TreeReader.hh
src/SignificanceAnalysis.o: ./include/base/TreeClassBase.h
src/SignificanceAnalysis.o: ./include/base/UserAnalysisBase.hh
src/SignificanceAnalysis.o: ./include/base/TreeReader.hh
src/SignificanceAnalysis.o: ./include/helper/pdgparticle.hh
src/SignificanceAnalysis.o: /usr/include/stdlib.h /usr/include/Availability.h
src/SignificanceAnalysis.o: /usr/include/AvailabilityInternal.h
src/SignificanceAnalysis.o: /usr/include/_types.h /usr/include/sys/_types.h
src/SignificanceAnalysis.o: /usr/include/sys/cdefs.h
src/SignificanceAnalysis.o: /usr/include/machine/_types.h
src/SignificanceAnalysis.o: /usr/include/i386/_types.h
src/SignificanceAnalysis.o: /usr/include/sys/wait.h /usr/include/sys/signal.h
src/SignificanceAnalysis.o: /usr/include/sys/appleapiopts.h
src/SignificanceAnalysis.o: /usr/include/machine/signal.h
src/SignificanceAnalysis.o: /usr/include/i386/signal.h
src/SignificanceAnalysis.o: /usr/include/i386/_structs.h
src/SignificanceAnalysis.o: /usr/include/sys/_structs.h
src/SignificanceAnalysis.o: /usr/include/machine/_structs.h
src/SignificanceAnalysis.o: /usr/include/sys/resource.h
src/SignificanceAnalysis.o: /usr/include/machine/endian.h
src/SignificanceAnalysis.o: /usr/include/i386/endian.h
src/SignificanceAnalysis.o: /usr/include/sys/_endian.h
src/SignificanceAnalysis.o: /usr/include/libkern/_OSByteOrder.h
src/SignificanceAnalysis.o: /usr/include/libkern/i386/_OSByteOrder.h
src/SignificanceAnalysis.o: /usr/include/alloca.h
src/SignificanceAnalysis.o: /usr/include/machine/types.h
src/SignificanceAnalysis.o: /usr/include/i386/types.h
src/SignificanceAnalysis.o: ./include/helper/Utilities.hh
src/SignificanceAnalysis.o: /usr/include/stdio.h /usr/include/secure/_stdio.h
src/SignificanceAnalysis.o: /usr/include/secure/_common.h
src/PhysQCAnalysis.o: ./include/base/TreeReader.hh
src/PhysQCAnalysis.o: ./include/base/TreeClassBase.h
src/PhysQCAnalysis.o: ./include/helper/AnaClass.hh /usr/include/math.h
src/PhysQCAnalysis.o: /usr/include/architecture/i386/math.h
src/PhysQCAnalysis.o: /usr/include/sys/cdefs.h /usr/include/stdlib.h
src/PhysQCAnalysis.o: /usr/include/Availability.h
src/PhysQCAnalysis.o: /usr/include/AvailabilityInternal.h
src/PhysQCAnalysis.o: /usr/include/_types.h /usr/include/sys/_types.h
src/PhysQCAnalysis.o: /usr/include/machine/_types.h
src/PhysQCAnalysis.o: /usr/include/i386/_types.h /usr/include/sys/wait.h
src/PhysQCAnalysis.o: /usr/include/sys/signal.h
src/PhysQCAnalysis.o: /usr/include/sys/appleapiopts.h
src/PhysQCAnalysis.o: /usr/include/machine/signal.h
src/PhysQCAnalysis.o: /usr/include/i386/signal.h /usr/include/i386/_structs.h
src/PhysQCAnalysis.o: /usr/include/sys/_structs.h
src/PhysQCAnalysis.o: /usr/include/machine/_structs.h
src/PhysQCAnalysis.o: /usr/include/sys/resource.h
src/PhysQCAnalysis.o: /usr/include/machine/endian.h
src/PhysQCAnalysis.o: /usr/include/i386/endian.h /usr/include/sys/_endian.h
src/PhysQCAnalysis.o: /usr/include/libkern/_OSByteOrder.h
src/PhysQCAnalysis.o: /usr/include/libkern/i386/_OSByteOrder.h
src/PhysQCAnalysis.o: /usr/include/alloca.h /usr/include/machine/types.h
src/PhysQCAnalysis.o: /usr/include/i386/types.h /usr/include/stdio.h
src/PhysQCAnalysis.o: /usr/include/secure/_stdio.h
src/PhysQCAnalysis.o: /usr/include/secure/_common.h /usr/include/time.h
src/PhysQCAnalysis.o: /usr/include/_structs.h ./include/helper/Utilities.hh
src/PhysQCAnalysis.o: ./include/helper/MetaTreeClassBase.h
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
src/RatioAnalysis.o: /usr/include/Availability.h
src/RatioAnalysis.o: /usr/include/AvailabilityInternal.h
src/RatioAnalysis.o: /usr/include/_types.h /usr/include/sys/_types.h
src/RatioAnalysis.o: /usr/include/sys/cdefs.h /usr/include/machine/_types.h
src/RatioAnalysis.o: /usr/include/i386/_types.h /usr/include/sys/wait.h
src/RatioAnalysis.o: /usr/include/sys/signal.h
src/RatioAnalysis.o: /usr/include/sys/appleapiopts.h
src/RatioAnalysis.o: /usr/include/machine/signal.h /usr/include/i386/signal.h
src/RatioAnalysis.o: /usr/include/i386/_structs.h /usr/include/sys/_structs.h
src/RatioAnalysis.o: /usr/include/machine/_structs.h
src/RatioAnalysis.o: /usr/include/sys/resource.h
src/RatioAnalysis.o: /usr/include/machine/endian.h /usr/include/i386/endian.h
src/RatioAnalysis.o: /usr/include/sys/_endian.h
src/RatioAnalysis.o: /usr/include/libkern/_OSByteOrder.h
src/RatioAnalysis.o: /usr/include/libkern/i386/_OSByteOrder.h
src/RatioAnalysis.o: /usr/include/alloca.h /usr/include/machine/types.h
src/RatioAnalysis.o: /usr/include/i386/types.h ./include/helper/Utilities.hh
src/RatioAnalysis.o: /usr/include/stdio.h /usr/include/secure/_stdio.h
src/RatioAnalysis.o: /usr/include/secure/_common.h
src/helper/TMctLib.o: ./include/helper/TMctLib.h ./include/helper/mctlib.h
src/helper/TMctLib.o: /usr/include/math.h
src/helper/TMctLib.o: /usr/include/architecture/i386/math.h
src/helper/TMctLib.o: /usr/include/sys/cdefs.h
src/helper/mctlib.o: ./include/helper/mctlib.h /usr/include/math.h
src/helper/mctlib.o: /usr/include/architecture/i386/math.h
src/helper/mctlib.o: /usr/include/sys/cdefs.h
src/helper/FPRatios.o: ./include/helper/FPRatios.hh
src/helper/AnaClass.o: ./include/helper/AnaClass.hh /usr/include/math.h
src/helper/AnaClass.o: /usr/include/architecture/i386/math.h
src/helper/AnaClass.o: /usr/include/sys/cdefs.h /usr/include/stdlib.h
src/helper/AnaClass.o: /usr/include/Availability.h
src/helper/AnaClass.o: /usr/include/AvailabilityInternal.h
src/helper/AnaClass.o: /usr/include/_types.h /usr/include/sys/_types.h
src/helper/AnaClass.o: /usr/include/machine/_types.h
src/helper/AnaClass.o: /usr/include/i386/_types.h /usr/include/sys/wait.h
src/helper/AnaClass.o: /usr/include/sys/signal.h
src/helper/AnaClass.o: /usr/include/sys/appleapiopts.h
src/helper/AnaClass.o: /usr/include/machine/signal.h
src/helper/AnaClass.o: /usr/include/i386/signal.h
src/helper/AnaClass.o: /usr/include/i386/_structs.h
src/helper/AnaClass.o: /usr/include/sys/_structs.h
src/helper/AnaClass.o: /usr/include/machine/_structs.h
src/helper/AnaClass.o: /usr/include/sys/resource.h
src/helper/AnaClass.o: /usr/include/machine/endian.h
src/helper/AnaClass.o: /usr/include/i386/endian.h /usr/include/sys/_endian.h
src/helper/AnaClass.o: /usr/include/libkern/_OSByteOrder.h
src/helper/AnaClass.o: /usr/include/libkern/i386/_OSByteOrder.h
src/helper/AnaClass.o: /usr/include/alloca.h /usr/include/machine/types.h
src/helper/AnaClass.o: /usr/include/i386/types.h /usr/include/stdio.h
src/helper/AnaClass.o: /usr/include/secure/_stdio.h
src/helper/AnaClass.o: /usr/include/secure/_common.h /usr/include/time.h
src/helper/AnaClass.o: /usr/include/_structs.h ./include/helper/Utilities.hh
src/helper/AnaClass.o: ./include/helper/MetaTreeClassBase.h
src/helper/AnaClass.o: ./include/helper/Utilities.hh
src/helper/Davismt2.o: ./include/helper/Davismt2.h /usr/include/math.h
src/helper/Davismt2.o: /usr/include/architecture/i386/math.h
src/helper/Davismt2.o: /usr/include/sys/cdefs.h
src/helper/LeptJetStat.o: ./include/helper/LeptJetStat.h
src/helper/Hemisphere.o: ./include/helper/Hemisphere.hh
src/helper/Hemisphere.o: ./include/helper/Utilities.hh /usr/include/stdio.h
src/helper/Hemisphere.o: /usr/include/_types.h /usr/include/sys/_types.h
src/helper/Hemisphere.o: /usr/include/sys/cdefs.h
src/helper/Hemisphere.o: /usr/include/machine/_types.h
src/helper/Hemisphere.o: /usr/include/i386/_types.h
src/helper/Hemisphere.o: /usr/include/secure/_stdio.h
src/helper/Hemisphere.o: /usr/include/secure/_common.h /usr/include/stdlib.h
src/helper/Hemisphere.o: /usr/include/Availability.h
src/helper/Hemisphere.o: /usr/include/AvailabilityInternal.h
src/helper/Hemisphere.o: /usr/include/sys/wait.h /usr/include/sys/signal.h
src/helper/Hemisphere.o: /usr/include/sys/appleapiopts.h
src/helper/Hemisphere.o: /usr/include/machine/signal.h
src/helper/Hemisphere.o: /usr/include/i386/signal.h
src/helper/Hemisphere.o: /usr/include/i386/_structs.h
src/helper/Hemisphere.o: /usr/include/sys/_structs.h
src/helper/Hemisphere.o: /usr/include/machine/_structs.h
src/helper/Hemisphere.o: /usr/include/sys/resource.h
src/helper/Hemisphere.o: /usr/include/machine/endian.h
src/helper/Hemisphere.o: /usr/include/i386/endian.h
src/helper/Hemisphere.o: /usr/include/sys/_endian.h
src/helper/Hemisphere.o: /usr/include/libkern/_OSByteOrder.h
src/helper/Hemisphere.o: /usr/include/libkern/i386/_OSByteOrder.h
src/helper/Hemisphere.o: /usr/include/alloca.h /usr/include/machine/types.h
src/helper/Hemisphere.o: /usr/include/i386/types.h
src/LeptJetMultAnalyzer.o: ./include/LeptJetMultAnalyzer.hh
src/LeptJetMultAnalyzer.o: ./include/base/TreeAnalyzerBase.hh
src/LeptJetMultAnalyzer.o: ./include/base/TreeReader.hh
src/LeptJetMultAnalyzer.o: ./include/base/TreeClassBase.h
src/LeptJetMultAnalyzer.o: ./include/helper/Utilities.hh /usr/include/stdio.h
src/LeptJetMultAnalyzer.o: /usr/include/_types.h /usr/include/sys/_types.h
src/LeptJetMultAnalyzer.o: /usr/include/sys/cdefs.h
src/LeptJetMultAnalyzer.o: /usr/include/machine/_types.h
src/LeptJetMultAnalyzer.o: /usr/include/i386/_types.h
src/LeptJetMultAnalyzer.o: /usr/include/secure/_stdio.h
src/LeptJetMultAnalyzer.o: /usr/include/secure/_common.h
src/LeptJetMultAnalyzer.o: /usr/include/stdlib.h /usr/include/Availability.h
src/LeptJetMultAnalyzer.o: /usr/include/AvailabilityInternal.h
src/LeptJetMultAnalyzer.o: /usr/include/sys/wait.h /usr/include/sys/signal.h
src/LeptJetMultAnalyzer.o: /usr/include/sys/appleapiopts.h
src/LeptJetMultAnalyzer.o: /usr/include/machine/signal.h
src/LeptJetMultAnalyzer.o: /usr/include/i386/signal.h
src/LeptJetMultAnalyzer.o: /usr/include/i386/_structs.h
src/LeptJetMultAnalyzer.o: /usr/include/sys/_structs.h
src/LeptJetMultAnalyzer.o: /usr/include/machine/_structs.h
src/LeptJetMultAnalyzer.o: /usr/include/sys/resource.h
src/LeptJetMultAnalyzer.o: /usr/include/machine/endian.h
src/LeptJetMultAnalyzer.o: /usr/include/i386/endian.h
src/LeptJetMultAnalyzer.o: /usr/include/sys/_endian.h
src/LeptJetMultAnalyzer.o: /usr/include/libkern/_OSByteOrder.h
src/LeptJetMultAnalyzer.o: /usr/include/libkern/i386/_OSByteOrder.h
src/LeptJetMultAnalyzer.o: /usr/include/alloca.h /usr/include/machine/types.h
src/LeptJetMultAnalyzer.o: /usr/include/i386/types.h
src/LeptJetMultAnalyzer.o: ./include/base/TreeReader.hh
src/LeptJetMultAnalyzer.o: ./include/MultiplicityAnalysis.hh
src/LeptJetMultAnalyzer.o: ./include/MultiplicityAnalysisBase.hh
src/LeptJetMultAnalyzer.o: ./include/base/UserAnalysisBase.hh
src/LeptJetMultAnalyzer.o: ./include/helper/pdgparticle.hh
src/LeptJetMultAnalyzer.o: ./include/helper/LeptJetStat.h
src/LeptJetMultAnalyzer.o: ./include/MassAnalysis.hh
src/LeptJetMultAnalyzer.o: ./include/helper/Davismt2.h /usr/include/math.h
src/LeptJetMultAnalyzer.o: /usr/include/architecture/i386/math.h
src/LeptJetMultAnalyzer.o: ./include/helper/TMctLib.h
src/LeptJetMultAnalyzer.o: ./include/helper/mctlib.h
src/LeptJetMultAnalyzer.o: ./include/helper/Hemisphere.hh
src/LeptJetMultAnalyzer.o: ./include/RatioAnalysis.hh
src/MuonPlotter.o: ./include/MuonPlotter.hh ./include/helper/AnaClass.hh
src/MuonPlotter.o: /usr/include/math.h /usr/include/architecture/i386/math.h
src/MuonPlotter.o: /usr/include/sys/cdefs.h /usr/include/stdlib.h
src/MuonPlotter.o: /usr/include/Availability.h
src/MuonPlotter.o: /usr/include/AvailabilityInternal.h /usr/include/_types.h
src/MuonPlotter.o: /usr/include/sys/_types.h /usr/include/machine/_types.h
src/MuonPlotter.o: /usr/include/i386/_types.h /usr/include/sys/wait.h
src/MuonPlotter.o: /usr/include/sys/signal.h /usr/include/sys/appleapiopts.h
src/MuonPlotter.o: /usr/include/machine/signal.h /usr/include/i386/signal.h
src/MuonPlotter.o: /usr/include/i386/_structs.h /usr/include/sys/_structs.h
src/MuonPlotter.o: /usr/include/machine/_structs.h
src/MuonPlotter.o: /usr/include/sys/resource.h /usr/include/machine/endian.h
src/MuonPlotter.o: /usr/include/i386/endian.h /usr/include/sys/_endian.h
src/MuonPlotter.o: /usr/include/libkern/_OSByteOrder.h
src/MuonPlotter.o: /usr/include/libkern/i386/_OSByteOrder.h
src/MuonPlotter.o: /usr/include/alloca.h /usr/include/machine/types.h
src/MuonPlotter.o: /usr/include/i386/types.h /usr/include/stdio.h
src/MuonPlotter.o: /usr/include/secure/_stdio.h /usr/include/secure/_common.h
src/MuonPlotter.o: /usr/include/time.h /usr/include/_structs.h
src/MuonPlotter.o: ./include/helper/Utilities.hh
src/MuonPlotter.o: ./include/helper/MetaTreeClassBase.h
src/MuonPlotter.o: ./include/helper/FPRatios.hh ./include/helper/Utilities.hh
src/MassPlotter.o: ./include/MassPlotter.hh
src/MassPlotter.o: ./include/helper/MassAnalysisTreeClass.h
src/MassPlotter.o: ./include/helper/Utilities.hh /usr/include/stdio.h
src/MassPlotter.o: /usr/include/_types.h /usr/include/sys/_types.h
src/MassPlotter.o: /usr/include/sys/cdefs.h /usr/include/machine/_types.h
src/MassPlotter.o: /usr/include/i386/_types.h /usr/include/secure/_stdio.h
src/MassPlotter.o: /usr/include/secure/_common.h /usr/include/stdlib.h
src/MassPlotter.o: /usr/include/Availability.h
src/MassPlotter.o: /usr/include/AvailabilityInternal.h
src/MassPlotter.o: /usr/include/sys/wait.h /usr/include/sys/signal.h
src/MassPlotter.o: /usr/include/sys/appleapiopts.h
src/MassPlotter.o: /usr/include/machine/signal.h /usr/include/i386/signal.h
src/MassPlotter.o: /usr/include/i386/_structs.h /usr/include/sys/_structs.h
src/MassPlotter.o: /usr/include/machine/_structs.h
src/MassPlotter.o: /usr/include/sys/resource.h /usr/include/machine/endian.h
src/MassPlotter.o: /usr/include/i386/endian.h /usr/include/sys/_endian.h
src/MassPlotter.o: /usr/include/libkern/_OSByteOrder.h
src/MassPlotter.o: /usr/include/libkern/i386/_OSByteOrder.h
src/MassPlotter.o: /usr/include/alloca.h /usr/include/machine/types.h
src/MassPlotter.o: /usr/include/i386/types.h /usr/include/math.h
src/MassPlotter.o: /usr/include/architecture/i386/math.h /usr/include/time.h
src/MassPlotter.o: /usr/include/_structs.h
src/helper/MetaTreeClassBase.o: ./include/helper/MetaTreeClassBase.h
src/helper/MassAnalysisTreeClass.o: ./include/helper/MassAnalysisTreeClass.h
src/JZBAnalyzer.o: ./include/JZBAnalyzer.hh
src/JZBAnalyzer.o: ./include/base/TreeAnalyzerBase.hh
src/JZBAnalyzer.o: ./include/base/TreeReader.hh
src/JZBAnalyzer.o: ./include/base/TreeClassBase.h
src/JZBAnalyzer.o: ./include/helper/Utilities.hh /usr/include/stdio.h
src/JZBAnalyzer.o: /usr/include/_types.h /usr/include/sys/_types.h
src/JZBAnalyzer.o: /usr/include/sys/cdefs.h /usr/include/machine/_types.h
src/JZBAnalyzer.o: /usr/include/i386/_types.h /usr/include/secure/_stdio.h
src/JZBAnalyzer.o: /usr/include/secure/_common.h /usr/include/stdlib.h
src/JZBAnalyzer.o: /usr/include/Availability.h
src/JZBAnalyzer.o: /usr/include/AvailabilityInternal.h
src/JZBAnalyzer.o: /usr/include/sys/wait.h /usr/include/sys/signal.h
src/JZBAnalyzer.o: /usr/include/sys/appleapiopts.h
src/JZBAnalyzer.o: /usr/include/machine/signal.h /usr/include/i386/signal.h
src/JZBAnalyzer.o: /usr/include/i386/_structs.h /usr/include/sys/_structs.h
src/JZBAnalyzer.o: /usr/include/machine/_structs.h
src/JZBAnalyzer.o: /usr/include/sys/resource.h /usr/include/machine/endian.h
src/JZBAnalyzer.o: /usr/include/i386/endian.h /usr/include/sys/_endian.h
src/JZBAnalyzer.o: /usr/include/libkern/_OSByteOrder.h
src/JZBAnalyzer.o: /usr/include/libkern/i386/_OSByteOrder.h
src/JZBAnalyzer.o: /usr/include/alloca.h /usr/include/machine/types.h
src/JZBAnalyzer.o: /usr/include/i386/types.h ./include/base/TreeReader.hh
src/JZBAnalyzer.o: ./include/JZBAnalysis.hh
src/JZBAnalyzer.o: ./include/base/UserAnalysisBase.hh
src/JZBAnalyzer.o: ./include/helper/pdgparticle.hh
src/JZBAnalyzer.o: ./include/helper/Monitor.hh
src/JZBAnalysis.o: ./include/helper/Utilities.hh /usr/include/stdio.h
src/JZBAnalysis.o: /usr/include/_types.h /usr/include/sys/_types.h
src/JZBAnalysis.o: /usr/include/sys/cdefs.h /usr/include/machine/_types.h
src/JZBAnalysis.o: /usr/include/i386/_types.h /usr/include/secure/_stdio.h
src/JZBAnalysis.o: /usr/include/secure/_common.h /usr/include/stdlib.h
src/JZBAnalysis.o: /usr/include/Availability.h
src/JZBAnalysis.o: /usr/include/AvailabilityInternal.h
src/JZBAnalysis.o: /usr/include/sys/wait.h /usr/include/sys/signal.h
src/JZBAnalysis.o: /usr/include/sys/appleapiopts.h
src/JZBAnalysis.o: /usr/include/machine/signal.h /usr/include/i386/signal.h
src/JZBAnalysis.o: /usr/include/i386/_structs.h /usr/include/sys/_structs.h
src/JZBAnalysis.o: /usr/include/machine/_structs.h
src/JZBAnalysis.o: /usr/include/sys/resource.h /usr/include/machine/endian.h
src/JZBAnalysis.o: /usr/include/i386/endian.h /usr/include/sys/_endian.h
src/JZBAnalysis.o: /usr/include/libkern/_OSByteOrder.h
src/JZBAnalysis.o: /usr/include/libkern/i386/_OSByteOrder.h
src/JZBAnalysis.o: /usr/include/alloca.h /usr/include/machine/types.h
src/JZBAnalysis.o: /usr/include/i386/types.h ./include/JZBAnalysis.hh
src/JZBAnalysis.o: ./include/base/TreeReader.hh
src/JZBAnalysis.o: ./include/base/TreeClassBase.h
src/JZBAnalysis.o: ./include/base/UserAnalysisBase.hh
src/JZBAnalysis.o: ./include/base/TreeReader.hh
src/JZBAnalysis.o: ./include/helper/pdgparticle.hh
src/JZBAnalysis.o: ./include/helper/Monitor.hh
src/SSDLAnalyzer.o: ./include/SSDLAnalyzer.hh
src/SSDLAnalyzer.o: ./include/base/TreeAnalyzerBase.hh
src/SSDLAnalyzer.o: ./include/base/TreeReader.hh
src/SSDLAnalyzer.o: ./include/base/TreeClassBase.h
src/SSDLAnalyzer.o: ./include/helper/Utilities.hh /usr/include/stdio.h
src/SSDLAnalyzer.o: /usr/include/_types.h /usr/include/sys/_types.h
src/SSDLAnalyzer.o: /usr/include/sys/cdefs.h /usr/include/machine/_types.h
src/SSDLAnalyzer.o: /usr/include/i386/_types.h /usr/include/secure/_stdio.h
src/SSDLAnalyzer.o: /usr/include/secure/_common.h /usr/include/stdlib.h
src/SSDLAnalyzer.o: /usr/include/Availability.h
src/SSDLAnalyzer.o: /usr/include/AvailabilityInternal.h
src/SSDLAnalyzer.o: /usr/include/sys/wait.h /usr/include/sys/signal.h
src/SSDLAnalyzer.o: /usr/include/sys/appleapiopts.h
src/SSDLAnalyzer.o: /usr/include/machine/signal.h /usr/include/i386/signal.h
src/SSDLAnalyzer.o: /usr/include/i386/_structs.h /usr/include/sys/_structs.h
src/SSDLAnalyzer.o: /usr/include/machine/_structs.h
src/SSDLAnalyzer.o: /usr/include/sys/resource.h /usr/include/machine/endian.h
src/SSDLAnalyzer.o: /usr/include/i386/endian.h /usr/include/sys/_endian.h
src/SSDLAnalyzer.o: /usr/include/libkern/_OSByteOrder.h
src/SSDLAnalyzer.o: /usr/include/libkern/i386/_OSByteOrder.h
src/SSDLAnalyzer.o: /usr/include/alloca.h /usr/include/machine/types.h
src/SSDLAnalyzer.o: /usr/include/i386/types.h ./include/helper/AnaClass.hh
src/SSDLAnalyzer.o: /usr/include/math.h /usr/include/architecture/i386/math.h
src/SSDLAnalyzer.o: /usr/include/time.h /usr/include/_structs.h
src/SSDLAnalyzer.o: ./include/helper/Utilities.hh
src/SSDLAnalyzer.o: ./include/helper/MetaTreeClassBase.h
src/SSDLAnalyzer.o: ./include/SSDLAnalysis.hh ./include/base/TreeReader.hh
src/SSDLAnalyzer.o: ./include/helper/Davismt2.h
src/SSDLAnalyzer.o: ./include/base/UserAnalysisBase.hh
src/SSDLAnalyzer.o: ./include/helper/pdgparticle.hh
src/SSDLAnalysis.o: ./include/SSDLAnalysis.hh /usr/include/time.h
src/SSDLAnalysis.o: /usr/include/_types.h /usr/include/sys/_types.h
src/SSDLAnalysis.o: /usr/include/sys/cdefs.h /usr/include/machine/_types.h
src/SSDLAnalysis.o: /usr/include/i386/_types.h /usr/include/_structs.h
src/SSDLAnalysis.o: /usr/include/sys/_structs.h
src/SSDLAnalysis.o: /usr/include/machine/_structs.h
src/SSDLAnalysis.o: /usr/include/i386/_structs.h
src/SSDLAnalysis.o: /usr/include/sys/appleapiopts.h
src/SSDLAnalysis.o: ./include/base/TreeReader.hh
src/SSDLAnalysis.o: ./include/base/TreeClassBase.h
src/SSDLAnalysis.o: ./include/helper/Davismt2.h /usr/include/math.h
src/SSDLAnalysis.o: /usr/include/architecture/i386/math.h
src/SSDLAnalysis.o: ./include/helper/Utilities.hh /usr/include/stdio.h
src/SSDLAnalysis.o: /usr/include/secure/_stdio.h
src/SSDLAnalysis.o: /usr/include/secure/_common.h /usr/include/stdlib.h
src/SSDLAnalysis.o: /usr/include/Availability.h
src/SSDLAnalysis.o: /usr/include/AvailabilityInternal.h
src/SSDLAnalysis.o: /usr/include/sys/wait.h /usr/include/sys/signal.h
src/SSDLAnalysis.o: /usr/include/machine/signal.h /usr/include/i386/signal.h
src/SSDLAnalysis.o: /usr/include/sys/resource.h /usr/include/machine/endian.h
src/SSDLAnalysis.o: /usr/include/i386/endian.h /usr/include/sys/_endian.h
src/SSDLAnalysis.o: /usr/include/libkern/_OSByteOrder.h
src/SSDLAnalysis.o: /usr/include/libkern/i386/_OSByteOrder.h
src/SSDLAnalysis.o: /usr/include/alloca.h /usr/include/machine/types.h
src/SSDLAnalysis.o: /usr/include/i386/types.h
src/SSDLAnalysis.o: ./include/base/UserAnalysisBase.hh
src/SSDLAnalysis.o: ./include/base/TreeReader.hh
src/SSDLAnalysis.o: ./include/helper/pdgparticle.hh
