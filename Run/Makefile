#
#   Makefile for complile of c++ code with call to ROOT libs
#
# CC -- the c compiler to use
# CXX -- the c++ compiler to use
# LD -- the linker to use
# CFLAGS -- compilation flag for c source files
# CXXFLAGS -- compilation flags for c++ source files
# CPPFLGS -- flags for the c-preprocessor (typically include file paths and symbols defined on the command line), used by c and c++
# LDFLAGS -- linker flags
# LDLIBS -- libraries to link


LDFLAGS	= 

GPP     = g++ -g -Wall

CC         = $(GPP)
CXX        = $(GPP)
LINKER	   = $(CC)

ROOTCFLAGS     = $(shell $(ROOTSYS)/bin/root-config --cflags) 
ROOTGLIBS      = $(shell $(ROOTSYS)/bin/root-config --glibs)  -lTreePlayer  -lMinuit -lMathMore -lRooFitCore -lRooFit -lXMLParser -lXMLIO -lCintex -lCint -lCore -lTMVA -lPyROOT
ROOTGLIBS     += -L$(ROOTCOREDIR)/lib -lTreePlayer  -lMinuit -lMathMore -lXMLParser -lXMLIO -lCintex -lCint -lCore -lTMVA -lPyROOT
ROOTGLIBS     += -L$(ROOTCOREDIR)/lib -lGoodRunsLists -lPileupReweighting -legammaAnalysisUtils -lMuonEfficiencyCorrections -lMuonMomentumCorrections -lHiggsZZ4lUtils -lTrigMuonEfficiency  -legammaFourMomentumError -lElectronEfficiencyCorrection -lPATCore -legParticleConditions -lZMassConstraint -lCLHEP -lApplyJetCalibration -lTileTripReader -lElectronPhotonSelectorTools -lElectronPhotonFourMomentumCorrection -lggFReweighting -lCategoriesMVA -legammaLayerRecalibTool -legammaMVACalib -lJHUReweighting -lH4lBDTWeights -lMadGraphME -lMyAnalysis -lD3PDReader 


INCDIRGR      = -IGoodRunsLists 
INCDIREA      = -IegammaAnalysisUtils
INCDIRME      = -IMuonEfficiencyCorrections
INCDIRMM      = -IMuonMomentumCorrections
INCDIRTM      = -ITrigMuonEfficiency
INCDIRPU      = -IPileupReweighting
INCDIRMI      = -IMuonIsolationCorrection
INCDIRHU      = -IHiggsZZ4lUtils
INCDIEGE      = -IegammaFourMomentumError
INCDIESF      = -IElectronEfficiencyCorrection
INCDIPAT      = -IPATCore
INCDIPCO      = -IegParticleConditions
INCDIRZM      = -IZMassConstraint
INCDIRCL      = -ICLHEP
INCDIRMA      = -IMyAnalysis
INCDIRDR      = -ID3PDReader
INCDIR        = -I$(ROOTCOREDIR)/include

CXXFLAGS       = -g $(ROOTCFLAGS) 
CXXFLAGS      += $(INCDIRGR)
CXXFLAGS      += $(INCDIRME)
CXXFLAGS      += $(INCDIRMM)
CXXFLAGS      += $(INCDIRTM)
CXXFLAGS      += $(INCDIREA)
CXXFLAGS      += $(INCDIRPU)
#CXXFLAGS      += $(INCDIRMI)
CXXFLAGS      += $(INCDIRHU)
CXXFLAGS      += $(INCDIEGE)
CXXFLAGS      += $(INCDIESF)
CXXFLAGS      += $(INCDIPAT)
CXXFLAGS      += $(INCDIPCO)
CXXFLAGS      += $(INCDIRZM)	
CXXFLAGS      += $(INCDIRCL)
CXXFLAGS      += $(INCDIRMA)
CXXFLAGS      += $(INCDIRDR)
CXXFLAGS      += $(INCDIR)
CXXFLAGS      += -DSTANDALONE # to pass the STANDALONE to the compiler
CXXFLAGS      += -D__USE_XOPEN2K8
CXXFLAGS      += -DROOTCORE
GLIBS         += $(ROOTGLIBS)


CCNAME = runAnalysis

OBJ  = $(CCNAME).o

SRC  = $(CCNAME).cxx

all:	runAnalysis


runAnalysis: runAnalysis.o
	$(LINKER) -o runAnalysis runAnalysis.o  $(GLIBS) $(SOFLAGS) 

runAnalysis.o: runAnalysis.cxx
	$(CC) -c runAnalysis.cxx $(CXXFLAGS) 

# suffix rule
.cxx.o:
	$(CXX) -c -pg $(CXXFLAGS) $(GDBFLAGS) $<

# clean
clean:
	rm -f *~ *.o *.o~ core



