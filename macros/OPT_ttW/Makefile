CC = g++
CFLAGS = -Wall -c -g

ROOFIT_INCLUDE := $(shell cd $(CMSSW_BASE); scram tool info roofitcore | grep INCLUDE= | sed 's|INCLUDE=||')
ROOFIT_LIBDIR := $(shell cd $(CMSSW_BASE); scram tool info roofitcore | grep LIBDIR= | sed 's|LIBDIR=||')

INCLUDES = -I. -I$(ROOTSYS)/include  -I$(ROOFIT_INCLUDE)/ -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools -I$(CMSSW_BASE)/src/UserCode/pandolf/ -I$(CMSSW_BASE)/src/DiLeptonAnalysis/NTupleProducer/macros

ROOTSYS  ?= ERROR_RootSysIsNotDefined

ROOTFLAG = `${ROOTSYS}/bin/root-config --cflags --libs`

EXTRALIBS  :=  -L$(ROOTSYS)/lib -L$(ROOFIT_LIBDIR)/ -lHtml -lMathCore -lGenVector -lMinuit -lEG -lRooFitCore -lRooFit -L$(CMSSW_RELEASE_BASE)/lib/$(SCRAM_ARCH) -lFWCoreFWLite -lCondFormatsJetMETObjects -lTreePlayer 



scanZBi: scanZBi.cpp DrawBase.o fitTools.o ZBiCalculator.o
	$(CC) -Wall $(INCLUDES) -o scanZBi scanZBi.cpp *.o ../src/SSDLPlotter.o ../src/SSDLDumper.o ../src/helper/*.o ../include/helper/BTagSFUtil/*.o $(ROOTFLAG) $(EXTRALIBS)

computeTTWZCrossSection: computeTTWZCrossSection.cpp DrawBase.o fitTools.o ZBiCalculator.o
	$(CC) -Wall $(INCLUDES) -o computeTTWZCrossSection computeTTWZCrossSection.cpp *.o ../src/*.o ../src/base/*.o ../src/helper/*.o ../include/helper/BTagSFUtil/*.o $(ROOTFLAG) $(EXTRALIBS)



DrawBase.o: $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/DrawBase.C fitTools.o
	$(CC) $(CFLAGS) $(INCLUDES) fitTools.o $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/DrawBase.C $(ROOTFLAG) $(EXTRALIBS) 

fitTools.o: $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/fitTools.C
	$(CC) $(CFLAGS) $(INCLUDES) $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/fitTools.C $(ROOTFLAG) $(EXTRALIBS)

ZBiCalculator.o: ZBiCalculator.cc
	$(CC) $(CFLAGS) $(INCLUDES) ZBiCalculator.cc $(ROOTFLAG) $(EXTRALIBS)

clean:
	rm *.o
