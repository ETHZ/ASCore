ROOTCFLAGS     = $(shell root-config --cflags)
ROOTLIBS       = $(shell root-config --libs)
ROOTGLIBS      = $(shell root-config --glibs)

CXX            = g++
CXXFLAGS       = -g -fPIC -Wno-deprecated -O -ansi -D_GNU_SOURCE -g -O2
LD             = g++
LDFLAGS        = -g
SOFLAGS        = -shared


ARCH         : = $(shell root-config --arch)
PLATFORM     : = $(shell root-config --platform)


CXXFLAGS      += $(ROOTCFLAGS)
#CXX           += -I./
LIBS           = $(ROOTLIBS) 

NGLIBS         = $(ROOTGLIBS) 
NGLIBS        += -lMinuit -lMinuit2
GLIBS          = $(filter-out -lNew, $(NGLIBS))

INCLUDEDIR1    = ./include/
CXX           += -I$(INCLUDEDIR1)
INCLUDEDIR2    = ./include/base/
CXX           += -I$(INCLUDEDIR2)
INCLUDEDIR3    = ./include/helper/
CXX           += -I$(INCLUDEDIR3)
OUTLIB         = ./obj/

.SUFFIXES: .cc,.C,.hh,.h
.PREFIXES: ./lib/

# Base Classes ===================================
TreeClassBase: src/base/TreeClassBase.C
	$(CXX) $(CXXFLAGS) -c -o $(OUTLIB)/TreeClassBase.o $<

TreeReader: src/base/TreeReader.cc
	$(CXX) $(CXXFLAGS) -c -o $(OUTLIB)/TreeReader.o $<

TreeAnalyzerBase: src/base/TreeAnalyzerBase.cc
	$(CXX) $(CXXFLAGS) -c -o $(OUTLIB)/TreeAnalyzerBase.o $<

UserAnalysisBase: src/base/UserAnalysisBase.cc
	$(CXX) $(CXXFLAGS) -c -o $(OUTLIB)/UserAnalysisBase.o $<

# Analyzers ======================================
TreeAnalyzer: src/TreeAnalyzer.cc
	$(CXX) $(CXXFLAGS) -c -o $(OUTLIB)/TreeAnalyzer.o $<

PhysQCAnalyzer: src/PhysQCAnalyzer.cc
	$(CXX) $(CXXFLAGS) -c -o $(OUTLIB)/PhysQCAnalyzer.o $<

# UserAnalysis Classes ===========================
DiLeptonAnalysis: src/DiLeptonAnalysis.cc
	$(CXX) $(CXXFLAGS) -c -o $(OUTLIB)/DiLeptonAnalysis.o $<

TreeCleaner: src/TreeCleaner.cc
	$(CXX) $(CXXFLAGS) -c -o $(OUTLIB)/TreeCleaner.o $<

MultiplicityAnalysis: src/MultiplicityAnalysis.cc
	$(CXX) $(CXXFLAGS) -c -o $(OUTLIB)/MultiplicityAnalysis.o $<

SignificanceAnalysis: src/SignificanceAnalysis.cc
	$(CXX) $(CXXFLAGS) -c -o $(OUTLIB)/SignificanceAnalysis.o $<

PhysQCAnalysis: src/PhysQCAnalysis.cc
	$(CXX) $(CXXFLAGS) -c -o $(OUTLIB)/PhysQCAnalysis.o $<

# Helper Classes =================================
AnaClass: src/helper/AnaClass.cc
	$(CXX) $(CXXFLAGS) -c -o $(OUTLIB)/AnaClass.o $<

Davismt2: src/helper/Davismt2.c
	$(CXX) $(CXXFLAGS) -c -o $(OUTLIB)/Davismt2.o $<

LeptJetStat: src/helper/LeptJetStat.cc
	$(CXX) $(CXXFLAGS) -c -o $(OUTLIB)/LeptJetStat.o $<

# Executables ====================================
RunTreeAnalyzer: src/exe/RunTreeAnalyzer.C
	$(CXX) $(CXXFLAGS) -ldl -o RunTreeAnalyzer $(OUTLIB)/*.o  $(GLIBS) $(LDFLAGS) $ $<

RunPhysQCAnalyzer: src/exe/RunPhysQCAnalyzer.C
	$(CXX) $(CXXFLAGS) -ldl -o RunPhysQCAnalyzer $(OUTLIB)/*.o  $(GLIBS) $(LDFLAGS) $ $<

clean:
	rm -f $(OUTLIB)*.o
	rm -f RunTreeAnalyzer
	rm -f RunPhysQCAnalyzer

PhysQC:
	touch src/exe/RunPhysQCAnalyzer.C
	make AnaClass
	make PhysQCAnalyzer
	make PhysQCAnalysis
	make RunPhysQCAnalyzer

all:
	touch src/exe/RunTreeAnalyzer.C
	touch src/exe/RunPhysQCAnalyzer.C
	make TreeClassBase
	make TreeReader
	make TreeAnalyzerBase
	make UserAnalysisBase
	make TreeAnalyzer
	make PhysQCAnalyzer
	make TreeCleaner
	make DiLeptonAnalysis
	make MultiplicityAnalysis
	make SignificanceAnalysis
	make PhysQCAnalysis
	make AnaClass
	make LeptJetStat
	make Davismt2
	make RunTreeAnalyzer
	make RunPhysQCAnalyzer
