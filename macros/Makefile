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
# NGLIBS        += -lMinuit -lMinuit2
NGLIBS        += -lMinuit
GLIBS          = $(filter-out -lNew, $(NGLIBS))

INCLUDEDIR     = ./include/
CXX           += -I$(INCLUDEDIR)
OUTLIB         = ./obj/
# OUTLIB         = $(INCLUDEDIR)

.SUFFIXES: .cc,.C,.hh,.h
.PREFIXES: ./lib/

#===================================
TreeClassBase: src/TreeClassBase.C
	$(CXX) $(CXXFLAGS) -c -o $(OUTLIB)/TreeClassBase.o $<

TreeReaderClass: src/TreeReader.cc
	$(CXX) $(CXXFLAGS) -c -o $(OUTLIB)/TreeReader.o $<

AnaClass: src/AnaClass.cc
	$(CXX) $(CXXFLAGS) -c -o $(OUTLIB)/AnaClass.o $<

Davismt2: src/Davismt2.c
	$(CXX) $(CXXFLAGS) -c -o $(OUTLIB)/Davismt2.o $<

# PtAveragePlots: src/PtAveragePlots.cc
# 	$(CXX) $(CXXFLAGS) -c -o $(OUTLIB)/PtAveragePlots.o $<

LeptJetStat: src/LeptJetStat.cc
	$(CXX) $(CXXFLAGS) -c -o $(OUTLIB)/LeptJetStat.o $<

RunTreeReader: src/RunTreeReader.C
	$(CXX) $(CXXFLAGS) -ldl -o RunTreeReader $(OUTLIB)/*.o  $(GLIBS) $(LDFLAGS) $ $<

# MakeAllPlots: src/MakeAllPlots.C
# 	$(CXX) $(CXXFLAGS) -ldl -o MakeAllPlots $(OUTLIB)/*.o  $(GLIBS) $(LDFLAGS) $ $<

# MakePlotList: src/MakePlotList.C
# 	$(CXX) $(CXXFLAGS) -ldl -o MakePlotList $(OUTLIB)/*.o  $(GLIBS) $(LDFLAGS) $ $<

# MakePtAvPlots: src/MakePtAvPlots.C
# 	$(CXX) $(CXXFLAGS) -ldl -o MakePtAvPlots $(OUTLIB)/*.o  $(GLIBS) $(LDFLAGS) $ $<

# MakePlotList2d: src/MakePlotList2d.C
# 	$(CXX) $(CXXFLAGS) -ldl -o MakePlotList2d $(OUTLIB)/*.o  $(GLIBS) $(LDFLAGS) $ $<

clean:
	rm -f $(OUTLIB)*.o
	rm -f RunTreeReader
	# rm -f MakeAllPlots
	# rm -f MakePlotList
	# rm -f MakePlotList2d

TreeReader:
	make AnaClass
	make TreeReaderClass
	touch src/RunTreeReader.C
	make RunTreeReader

all:
	make TreeClassBase
	make TreeReaderClass
	make AnaClass
	make LeptJetStat
	make Davismt2
	# make PtAveragePlots
	make RunTreeReader
	# make MakeAllPlots
	# make MakePlotList
	# make MakePtAvPlots
	# make MakePlotList2d
