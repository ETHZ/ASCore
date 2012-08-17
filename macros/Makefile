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

NGLIBS         = $(ROOTGLIBS) -lMinuit -lMinuit2 -lTreePlayer -lGeom
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
	$(RM) RunDiPhotonJetsAnalyzer

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
src/base/TreeAnalyzerBase.o: /usr/include/sys/_symbol_aliasing.h
src/base/TreeAnalyzerBase.o: /usr/include/sys/_posix_availability.h
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
src/base/TreeAnalyzerBase.o: /usr/include/assert.h
src/base/UserAnalysisBase.o: ./include/base/TreeReader.hh
src/base/UserAnalysisBase.o: ./include/base/TreeClassBase.h
src/base/UserAnalysisBase.o: /usr/include/stdlib.h
src/base/UserAnalysisBase.o: /usr/include/Availability.h
src/base/UserAnalysisBase.o: /usr/include/AvailabilityInternal.h
src/base/UserAnalysisBase.o: /usr/include/_types.h /usr/include/sys/_types.h
src/base/UserAnalysisBase.o: /usr/include/sys/cdefs.h
src/base/UserAnalysisBase.o: /usr/include/sys/_symbol_aliasing.h
src/base/UserAnalysisBase.o: /usr/include/sys/_posix_availability.h
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
src/base/UserAnalysisBase.o: ./include/helper/Monitor.hh
src/base/UserAnalysisBase.o: ./include/base/UserAnalysisBase.hh
src/base/UserAnalysisBase.o: ./include/base/TreeReader.hh
src/base/UserAnalysisBase.o: ./include/helper/Utilities.hh
src/base/UserAnalysisBase.o: /usr/include/stdio.h
src/base/UserAnalysisBase.o: /usr/include/secure/_stdio.h
src/base/UserAnalysisBase.o: /usr/include/secure/_common.h
src/base/UserAnalysisBase.o: ./include/helper/PUWeight.h
src/base/UserAnalysisBase.o: ./include/helper/Lumi3DReWeighting_standalone.hh
src/UserAnalyzer.o: ./include/UserAnalyzer.hh
src/UserAnalyzer.o: ./include/base/TreeAnalyzerBase.hh
src/UserAnalyzer.o: ./include/base/TreeReader.hh
src/UserAnalyzer.o: ./include/base/TreeClassBase.h
src/UserAnalyzer.o: ./include/helper/Utilities.hh /usr/include/stdio.h
src/UserAnalyzer.o: /usr/include/sys/cdefs.h
src/UserAnalyzer.o: /usr/include/sys/_symbol_aliasing.h
src/UserAnalyzer.o: /usr/include/sys/_posix_availability.h
src/UserAnalyzer.o: /usr/include/Availability.h
src/UserAnalyzer.o: /usr/include/AvailabilityInternal.h /usr/include/_types.h
src/UserAnalyzer.o: /usr/include/sys/_types.h /usr/include/machine/_types.h
src/UserAnalyzer.o: /usr/include/i386/_types.h /usr/include/secure/_stdio.h
src/UserAnalyzer.o: /usr/include/secure/_common.h /usr/include/stdlib.h
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
src/UserAnalyzer.o: ./include/helper/PUWeight.h
src/UserAnalyzer.o: ./include/helper/Lumi3DReWeighting_standalone.hh
src/TreeSkimmer.o: ./include/TreeSkimmer.hh
src/TreeSkimmer.o: ./include/base/TreeAnalyzerBase.hh
src/TreeSkimmer.o: ./include/base/TreeReader.hh
src/TreeSkimmer.o: ./include/base/TreeClassBase.h
src/TreeSkimmer.o: ./include/helper/Utilities.hh /usr/include/stdio.h
src/TreeSkimmer.o: /usr/include/sys/cdefs.h
src/TreeSkimmer.o: /usr/include/sys/_symbol_aliasing.h
src/TreeSkimmer.o: /usr/include/sys/_posix_availability.h
src/TreeSkimmer.o: /usr/include/Availability.h
src/TreeSkimmer.o: /usr/include/AvailabilityInternal.h /usr/include/_types.h
src/TreeSkimmer.o: /usr/include/sys/_types.h /usr/include/machine/_types.h
src/TreeSkimmer.o: /usr/include/i386/_types.h /usr/include/secure/_stdio.h
src/TreeSkimmer.o: /usr/include/secure/_common.h /usr/include/stdlib.h
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
src/UserAnalysis.o: ./include/helper/Utilities.hh /usr/include/stdio.h
src/UserAnalysis.o: /usr/include/sys/cdefs.h
src/UserAnalysis.o: /usr/include/sys/_symbol_aliasing.h
src/UserAnalysis.o: /usr/include/sys/_posix_availability.h
src/UserAnalysis.o: /usr/include/Availability.h
src/UserAnalysis.o: /usr/include/AvailabilityInternal.h /usr/include/_types.h
src/UserAnalysis.o: /usr/include/sys/_types.h /usr/include/machine/_types.h
src/UserAnalysis.o: /usr/include/i386/_types.h /usr/include/secure/_stdio.h
src/UserAnalysis.o: /usr/include/secure/_common.h /usr/include/stdlib.h
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
src/UserAnalysis.o: ./include/helper/PUWeight.h
src/UserAnalysis.o: ./include/helper/Lumi3DReWeighting_standalone.hh
src/helper/PUWeight.o: ./include/helper/PUWeight.h
src/helper/AnaClass.o: ./include/helper/AnaClass.hh /usr/include/math.h
src/helper/AnaClass.o: /usr/include/architecture/i386/math.h
src/helper/AnaClass.o: /usr/include/sys/cdefs.h
src/helper/AnaClass.o: /usr/include/sys/_symbol_aliasing.h
src/helper/AnaClass.o: /usr/include/sys/_posix_availability.h
src/helper/AnaClass.o: /usr/include/stdlib.h /usr/include/Availability.h
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
src/helper/Davismt2.o: /usr/include/sys/_symbol_aliasing.h
src/helper/Davismt2.o: /usr/include/sys/_posix_availability.h
src/helper/LeptJetStat.o: ./include/helper/LeptJetStat.h
src/helper/Hemisphere.o: ./include/helper/Hemisphere.hh
src/helper/Hemisphere.o: ./include/helper/Utilities.hh /usr/include/stdio.h
src/helper/Hemisphere.o: /usr/include/sys/cdefs.h
src/helper/Hemisphere.o: /usr/include/sys/_symbol_aliasing.h
src/helper/Hemisphere.o: /usr/include/sys/_posix_availability.h
src/helper/Hemisphere.o: /usr/include/Availability.h
src/helper/Hemisphere.o: /usr/include/AvailabilityInternal.h
src/helper/Hemisphere.o: /usr/include/_types.h /usr/include/sys/_types.h
src/helper/Hemisphere.o: /usr/include/machine/_types.h
src/helper/Hemisphere.o: /usr/include/i386/_types.h
src/helper/Hemisphere.o: /usr/include/secure/_stdio.h
src/helper/Hemisphere.o: /usr/include/secure/_common.h /usr/include/stdlib.h
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
src/helper/MetaTreeClassBase.o: ./include/helper/MetaTreeClassBase.h
src/helper/Lumi3DReWeighting_standalone.o: ./include/helper/Lumi3DReWeighting_standalone.hh
src/EnergyCorrection.o: ./include/EnergyCorrection.hh
src/EnergyCorrection.o: ./include/base/TreeReader.hh
src/EnergyCorrection.o: ./include/base/TreeClassBase.h /usr/include/assert.h
src/EnergyCorrection.o: /usr/include/sys/cdefs.h
src/EnergyCorrection.o: /usr/include/sys/_symbol_aliasing.h
src/EnergyCorrection.o: /usr/include/sys/_posix_availability.h
src/DiPhotonMiniTree.o: ./include/helper/Utilities.hh /usr/include/stdio.h
src/DiPhotonMiniTree.o: /usr/include/sys/cdefs.h
src/DiPhotonMiniTree.o: /usr/include/sys/_symbol_aliasing.h
src/DiPhotonMiniTree.o: /usr/include/sys/_posix_availability.h
src/DiPhotonMiniTree.o: /usr/include/Availability.h
src/DiPhotonMiniTree.o: /usr/include/AvailabilityInternal.h
src/DiPhotonMiniTree.o: /usr/include/_types.h /usr/include/sys/_types.h
src/DiPhotonMiniTree.o: /usr/include/machine/_types.h
src/DiPhotonMiniTree.o: /usr/include/i386/_types.h
src/DiPhotonMiniTree.o: /usr/include/secure/_stdio.h
src/DiPhotonMiniTree.o: /usr/include/secure/_common.h /usr/include/stdlib.h
src/DiPhotonMiniTree.o: /usr/include/sys/wait.h /usr/include/sys/signal.h
src/DiPhotonMiniTree.o: /usr/include/sys/appleapiopts.h
src/DiPhotonMiniTree.o: /usr/include/machine/signal.h
src/DiPhotonMiniTree.o: /usr/include/i386/signal.h
src/DiPhotonMiniTree.o: /usr/include/i386/_structs.h
src/DiPhotonMiniTree.o: /usr/include/sys/_structs.h
src/DiPhotonMiniTree.o: /usr/include/machine/_structs.h
src/DiPhotonMiniTree.o: /usr/include/sys/resource.h
src/DiPhotonMiniTree.o: /usr/include/machine/endian.h
src/DiPhotonMiniTree.o: /usr/include/i386/endian.h /usr/include/sys/_endian.h
src/DiPhotonMiniTree.o: /usr/include/libkern/_OSByteOrder.h
src/DiPhotonMiniTree.o: /usr/include/libkern/i386/_OSByteOrder.h
src/DiPhotonMiniTree.o: /usr/include/alloca.h /usr/include/machine/types.h
src/DiPhotonMiniTree.o: /usr/include/i386/types.h
src/DiPhotonMiniTree.o: ./include/DiPhotonMiniTree.hh
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
src/DiPhotonPurity.o: /usr/include/sys/cdefs.h
src/DiPhotonPurity.o: /usr/include/sys/_symbol_aliasing.h
src/DiPhotonPurity.o: /usr/include/sys/_posix_availability.h
src/DiPhotonPurity.o: /usr/include/Availability.h
src/DiPhotonPurity.o: /usr/include/AvailabilityInternal.h
src/DiPhotonPurity.o: /usr/include/_types.h /usr/include/sys/_types.h
src/DiPhotonPurity.o: /usr/include/machine/_types.h
src/DiPhotonPurity.o: /usr/include/i386/_types.h /usr/include/secure/_stdio.h
src/DiPhotonPurity.o: /usr/include/secure/_common.h /usr/include/stdlib.h
src/DiPhotonPurity.o: /usr/include/sys/wait.h /usr/include/sys/signal.h
src/DiPhotonPurity.o: /usr/include/sys/appleapiopts.h
src/DiPhotonPurity.o: /usr/include/machine/signal.h
src/DiPhotonPurity.o: /usr/include/i386/signal.h /usr/include/i386/_structs.h
src/DiPhotonPurity.o: /usr/include/sys/_structs.h
src/DiPhotonPurity.o: /usr/include/machine/_structs.h
src/DiPhotonPurity.o: /usr/include/sys/resource.h
src/DiPhotonPurity.o: /usr/include/machine/endian.h
src/DiPhotonPurity.o: /usr/include/i386/endian.h /usr/include/sys/_endian.h
src/DiPhotonPurity.o: /usr/include/libkern/_OSByteOrder.h
src/DiPhotonPurity.o: /usr/include/libkern/i386/_OSByteOrder.h
src/DiPhotonPurity.o: /usr/include/alloca.h /usr/include/machine/types.h
src/DiPhotonPurity.o: /usr/include/i386/types.h ./include/DiPhotonPurity.hh
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
src/DiPhotonJetsAnalyzer.o: /usr/include/stdio.h /usr/include/sys/cdefs.h
src/DiPhotonJetsAnalyzer.o: /usr/include/sys/_symbol_aliasing.h
src/DiPhotonJetsAnalyzer.o: /usr/include/sys/_posix_availability.h
src/DiPhotonJetsAnalyzer.o: /usr/include/Availability.h
src/DiPhotonJetsAnalyzer.o: /usr/include/AvailabilityInternal.h
src/DiPhotonJetsAnalyzer.o: /usr/include/_types.h /usr/include/sys/_types.h
src/DiPhotonJetsAnalyzer.o: /usr/include/machine/_types.h
src/DiPhotonJetsAnalyzer.o: /usr/include/i386/_types.h
src/DiPhotonJetsAnalyzer.o: /usr/include/secure/_stdio.h
src/DiPhotonJetsAnalyzer.o: /usr/include/secure/_common.h
src/DiPhotonJetsAnalyzer.o: /usr/include/stdlib.h /usr/include/sys/wait.h
src/DiPhotonJetsAnalyzer.o: /usr/include/sys/signal.h
src/DiPhotonJetsAnalyzer.o: /usr/include/sys/appleapiopts.h
src/DiPhotonJetsAnalyzer.o: /usr/include/machine/signal.h
src/DiPhotonJetsAnalyzer.o: /usr/include/i386/signal.h
src/DiPhotonJetsAnalyzer.o: /usr/include/i386/_structs.h
src/DiPhotonJetsAnalyzer.o: /usr/include/sys/_structs.h
src/DiPhotonJetsAnalyzer.o: /usr/include/machine/_structs.h
src/DiPhotonJetsAnalyzer.o: /usr/include/sys/resource.h
src/DiPhotonJetsAnalyzer.o: /usr/include/machine/endian.h
src/DiPhotonJetsAnalyzer.o: /usr/include/i386/endian.h
src/DiPhotonJetsAnalyzer.o: /usr/include/sys/_endian.h
src/DiPhotonJetsAnalyzer.o: /usr/include/libkern/_OSByteOrder.h
src/DiPhotonJetsAnalyzer.o: /usr/include/libkern/i386/_OSByteOrder.h
src/DiPhotonJetsAnalyzer.o: /usr/include/alloca.h
src/DiPhotonJetsAnalyzer.o: /usr/include/machine/types.h
src/DiPhotonJetsAnalyzer.o: /usr/include/i386/types.h
src/DiPhotonJetsAnalyzer.o: ./include/base/TreeReader.hh
src/DiPhotonJetsAnalyzer.o: ./include/DiPhotonPurity.hh
src/DiPhotonJetsAnalyzer.o: ./include/base/UserAnalysisBase.hh
src/DiPhotonJetsAnalyzer.o: ./include/helper/pdgparticle.hh
src/DiPhotonJetsAnalyzer.o: ./include/helper/PUWeight.h
src/DiPhotonJetsAnalyzer.o: ./include/helper/Lumi3DReWeighting_standalone.hh
src/DiPhotonJetsAnalyzer.o: ./include/EnergyCorrection.hh
src/DiPhotonJetsAnalyzer.o: ./include/DiPhotonMiniTree.hh
