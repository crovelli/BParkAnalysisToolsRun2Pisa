ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs) -lTMVA

CXX           = g++
CXXFLAGS      = -g -fPIC -Wno-deprecated -O -ansi -D_GNU_SOURCE -g -O2
LD            = g++
LDFLAGS       = -g
SOFLAGS       = -shared


ARCH         := $(shell root-config --arch)
PLATFORM     := $(shell root-config --platform)


CXXFLAGS      += $(ROOTCFLAGS)
#CXX           += -I./
LIBS           = $(ROOTLIBS) 

NGLIBS         = $(ROOTGLIBS) 
GLIBS          = $(filter-out -lNew, $(NGLIBS))

INCLUDEDIR       = ./
CXX	         += -I$(INCLUDEDIR) -I.
OUTLIB	         = $(INCLUDEDIR)/lib/

.SUFFIXES: .cc,.C,.hh,.h
.PREFIXES: ./lib/


$(OUTLIB)BParkBase.o: $(INCLUDEDIR)/src/BParkBase.C \
	$(INCLUDEDIR)/src/BPark.cc \
	$(INCLUDEDIR)/src/SkimmerPisa.cc 
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)BParkBase.o $<
$(OUTLIB)BPark.o: $(INCLUDEDIR)/src/BPark.cc 
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)BPark.o $<
$(OUTLIB)SkimmerPisa.o: $(INCLUDEDIR)/src/SkimmerPisa.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)SkimmerPisa.o $<

# ==================== BParkApp =============================================
BParkApp:  $(INCLUDEDIR)/src/BParkApp.C \
	$(OUTLIB)BParkBase.o \
	$(OUTLIB)BPark.o 
	$(CXX) $(CXXFLAGS) -ldl -o BParkApp $(OUTLIB)/*.o $(GLIBS) $(LDFLAGS) $ $<
BParkApp.clean:
	rm -f BParkApp

# ==================== reduced trees =============================================

clean:
	rm -f $(OUTLIB)*.o
	rm -f BParkApp

all:  BParkApp
