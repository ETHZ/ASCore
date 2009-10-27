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

INCLUDEDIR     = ./include/
CXX           += -I$(INCLUDEDIR)
OUTLIB         = ./obj/

.SUFFIXES: .cc,.C,.hh,.h
.PREFIXES: ./lib/

#===================================
TreeClassBase: src/TreeClassBase.C
	$(CXX) $(CXXFLAGS) -c -o $(OUTLIB)/TreeClassBase.o $<

TreeReader: src/TreeReader.cc
	$(CXX) $(CXXFLAGS) -c -o $(OUTLIB)/TreeReader.o $<

AnaClass: src/AnaClass.cc
	$(CXX) $(CXXFLAGS) -c -o $(OUTLIB)/AnaClass.o $<

Davismt2: src/Davismt2.c
	$(CXX) $(CXXFLAGS) -c -o $(OUTLIB)/Davismt2.o $<

LeptJetStat: src/LeptJetStat.cc
	$(CXX) $(CXXFLAGS) -c -o $(OUTLIB)/LeptJetStat.o $<

RunTreeReader: src/RunTreeReader.C
	$(CXX) $(CXXFLAGS) -ldl -o RunTreeReader $(OUTLIB)/*.o  $(GLIBS) $(LDFLAGS) $ $<

clean:
	rm -f $(OUTLIB)*.o
	rm -f RunTreeReader

Reader:
	make AnaClass
	make TreeReader
	touch src/RunTreeReader.C
	make RunTreeReader

all:
	touch src/RunTreeReader.C
	make TreeClassBase
	make TreeReader
	make AnaClass
	make LeptJetStat
	make Davismt2
	make RunTreeReader
