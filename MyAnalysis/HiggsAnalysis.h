#ifndef HIGGANALYSIS_H
#define HIGGANALYSIS_H

#include "D3PDReader/Event.h"
#include "D3PDReader/MuonD3PDObject.h"
#include "D3PDReader/ElectronD3PDObject.h"
#include "D3PDReader/TruthParticleD3PDObject.h"

#include "MyAnalysis/macroDef.h"
#include "MyAnalysis/ChargedLepton.h"
#include "MyAnalysis/MuonObject.h"
#include "MyAnalysis/ElectronObject.h"
#include "MyAnalysis/JetsObject.h"
#include "MyAnalysis/DiLepton.h"
#include "MyAnalysis/QuadLepton.h"
#include "MyAnalysis/HistContainer.h"
#include "MyAnalysis/Correction.h"
#include "MyAnalysis/OutputTree.h"
#include "MyAnalysis/OutputTreeSys.h"

#include "PileupReweighting/TPileupReweighting.h"
#include "egammaAnalysisUtils/CaloIsoCorrection.h"
#include "egammaAnalysisUtils/VertexPositionReweightingTool.h"

#include "Math/Interpolator.h"
#include "Math/InterpolationTypes.h"
#include "Math/Polynomial.h"

#include "TrigMuonEfficiency/MuonTriggerMatching.h"
#include "TrigMuonEfficiency/ElectronTriggerMatching.h"
#include "TrigMuonEfficiency/TriggerNavigationVariables.h"
#include "TrigMuonEfficiency/LeptonTriggerSF.h"

#include "egammaAnalysisUtils/FsrPhotons.h"
#include "egammaFourMomentumError/egammaFourMomentumError.h" 
#include "ZMassConstraint/ConstraintFit.h"
#include "ZMassConstraint/ConstraintFitInput.h"
#include "ZMassConstraint/ConstraintFitOutput.h"
#include "ZMassConstraint/CovMatrixTools.h"

#include "HiggsZZ4lUtils/McOverlapRemoval.h"
#include "HiggsZZ4lUtils/GetElicityAngles.h"
#include "HiggsZZ4lUtils/BRCorrection.h"

#include "GoodRunsLists/TGoodRunsListReader.h"
#include "GoodRunsLists/TGoodRunsListWriter.h"
#include "GoodRunsLists/TGRLCollection.h"
#include "GoodRunsLists/TGoodRunsList.h"

#include "TileTripReader/TTileTripReader.h"

#include "HiggsZZ4lUtils/HiggsCrossSection.h"
#include "HiggsZZ4lUtils/BkgCrossSection.h"
#include "HiggsZZ4lUtils/H4lBrRatio.h"

#include "ggFReweighting/ggFReweighting.h"

#include "CategoriesMVA/CategoriesMVA.h"

#include "JHUReweighting/JHUPtReweighting.h"

#include "H4lBDTWeights/H4lBDTWeights.h"

#include "BCHCleaningTool/BCHCleaningToolRoot.h"

#include <TTree.h>
#include <TChain.h>
#include <TString.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TLorentzVector.h>

#include <vector>
#include <algorithm>
#include "math.h"


using namespace std;
class HiggsAnalysis 
{
	public :
		// Main Variables
		D3PDReader::Event *event;
		TTree *physicsTree; 

		MuonObject *muonObj;
		ElectronObject *electronObj;
		ElectronObject *electronObjLoose;
		JetsObject *jetsObj;
		JetsObject *jetsObjTruth;
		JetsObject *jetsObj_Fid;
		JetsObject *jetsObjTruth_Fid;

		// Final Higgs Quadleptong
		QuadLepton *higgsCandidate4Mu;
		QuadLepton *higgsCandidate4El;
		QuadLepton *higgsCandidate2L2L;

		// Helpers to keep track of what analysis to do
		Bool_t do4Mu;
		Bool_t do4El;
		Bool_t do2L2L;
		Bool_t dollee;
		Bool_t doCR;
		Bool_t doSysTree;
	
		// For Moriond CutFlow
		Bool_t doMoriondConfig;

		// To Control all the smearing
		Bool_t doCorr;
		Bool_t doSmearD0Z0;
		Bool_t doSmearMC;
		Bool_t doScaleMuon;		
		Bool_t doScaleData;
		Bool_t doScaleCrack;
		Bool_t doCaloIsolationCorr;
		Bool_t doJetCal;
		Bool_t doEPcombination;
		
		// To Control all the weight
		Bool_t doWeight;		
		Bool_t doPileupWeight;
		Bool_t doZVertexWeight;
		Bool_t doggFWeight;		
		Bool_t doJHUWeight;		
		Bool_t doTTBarOverlap;
		Bool_t doZZOverlap;
		Bool_t doZBBOverlap;
		Bool_t doScaleEfficiency;
		Bool_t doTriggerEfficiency;

		// To Control if I want to use the new Likelihood Tool
		Bool_t useLikelihood;

		// To Kepp track it is it a debug call
		Bool_t isDebugCall;
		Bool_t isDebugCallllee;
		Bool_t printWeight;

		// BCHCut
		Bool_t temp_BCHCutMedium;
		Bool_t temp_BCHCutTight;
		
		// Tools
		Root::TPileupReweighting* pileupTool;
  		
		TriggerNavigationVariables* triggerNavigationVariables;
  		MuonTriggerMatching *muonTriggerMatchTool;
  		ElectronTriggerMatching *electronTriggerMatchTool;
 		LeptonTriggerSF *leptonSF;
 		
		VertexPositionReweightingTool *zVertexTool;
		
		Root::TGoodRunsList grl;
		
		Root::TTileTripReader* TileTrip;
		BCHTool::BCHCleaningToolRoot* thebchTool;
		BCHTool::BCHCleaningToolRoot* thebchToolMedium;
		
		HiggsCrossSection Higgs_xs;
		H4lBrRatio *higgs_bror;
		BRCorrection * brCorr;

		ggFReweighting *ggFReweight;

  		JHUPtReweighting *ptJHUReweighting;
		
		//// BDT for the categories
  		CategoriesMVA* CategoriesDiscriminantTool;

		// BDT for 2D fit
		H4lBDTWeights *BDTtool;

		// All the corrections
		Correction *corr;

		// Histrogram Container
		HistContainer *Hist;

		// For Output Tree
		OutputTree *outputTree;
		OutputTree *outputTreeCR;
		OutputTree *outputTreelleeCR;
		OutputTreeSys *outputTreeSys;
		TString outputFilePath;

		// Helper Variables for event
		Int_t curEvent;
		Int_t dataYear;
		Bool_t is2012;
		Int_t streamName;
		TString dataPeriod;
		Bool_t noTauSample;
		Int_t generatorName;
		Int_t curMCCollection;
		Int_t curDataCollection;
		Int_t curCalibration;
		Int_t anaType;

		Int_t electronCollection;
		Int_t muonCollection;
		Int_t sampleProdType;

		Double_t CM_E;
		Bool_t isMC;

		Int_t runNumber_sf;
		Int_t lbn_sf;
		Int_t mcRunNumber;

		TString currFileName;
		TString currSampleName;
		TH1F * countingHist;

		Bool_t printProductionTag;

		// For running on the gird
		Bool_t runningGrid;
		TString gridFileName;

		// To choose the calibration for the data
		Bool_t useNewGeoData;

		// Variables for event list
		Bool_t printEventList; // need to set this up for it to actuall print
		Bool_t printMass; // If I want to print the mass for each event
		ofstream file4Mu;
		ofstream file4El;
		ofstream file2El2Mu;
		ofstream file2Mu2El;
		ofstream file2L2L;
		ofstream fileSampleName;

		// For Trigger Matching
		Bool_t D3DPTriggerDev;
		
		// To keep track of the electron set for 2011 and 2012
		D3PDReader::ElectronD3PDObject * el_cur;

		// Varibles that store the final physics objects
		vector<ChargedLepton *> muEvent;
		vector<ChargedLepton *> elEvent;
		vector<ChargedLepton *> elLooseEvent;
		vector<ChargedLepton *> jetsEvent;
		vector<ChargedLepton *> jetsTruthEvent;
		vector<ChargedLepton *> jetsEvent_Fid;
		vector<ChargedLepton *> jetsTruthEvent_Fid;
	
		// Physics objects before final overlap removal
		vector<ChargedLepton *> muOverlap;
		vector<ChargedLepton *> elOverlap;
		vector<ChargedLepton *> elLooseOverlap;
		vector<ChargedLepton *> jetsOverlap;
		vector<ChargedLepton *> jetsOverlap_Fid;

		vector<DiLepton *> diEvent2Mu;
		vector<DiLepton *> diEvent2El;
		vector<QuadLepton *> eventToClean4mu;
		vector<QuadLepton *> eventToClean4el;
		vector<QuadLepton *> eventToClean2mu2el;
		vector<QuadLepton *> eventToCleanCR;
		vector<DiLepton *> eventToCleanDi;

		// Varibles to Count number of events that have passed
		// Plus Histrograms to fill
		Int_t nCut;
		TString *cutName;
		Int_t *cutPass;
		Double_t *cutPassW;

		Int_t nMuCut;
		TString *cutMuName;
		Int_t *cutMuPass;

		Int_t nElCut;
		TString *cutElName;
		Int_t *cutElPass;
		Int_t *cutElLoosePass;

		Int_t nJetsCut;
		TString *cutJetsName;
		Int_t *cutJetsPass;

		Int_t nCH;
		TString *cutCHName;
		
		Int_t *cut4MuPass;
		Int_t *cut4ElPass;
		Int_t *cut2L2LPass;
		Int_t *cutlleePass;
		Double_t *cut4MuPassW;
		Double_t *cut4ElPassW;
		Double_t *cut2L2LPassW;

		Int_t nListTruthQuadType;
		TString *truthQuadType;
		Int_t *nTruthQuadType;

		Int_t nProdCH;
		TString *prodCHName;
		Int_t *prodCH4Mu;
		Int_t *prodCH4El;
		Int_t *prodCH2L2L;
		Int_t *prodCHllee;

		// For Categories
		Float_t dijet_invmass;
		Float_t dijet_deltaeta;
		Float_t leading_jet_pt;
		Float_t leading_jet_eta;
		Float_t subleading_jet_pt;

		Float_t BDT_discriminant_VBF;
		Float_t BDT_discriminant_HadVH;
		
		// Constructors and Destructors
		HiggsAnalysis(TTree *physics, Bool_t tRunningGrid = false, TString tgridFileName = "runningGridMismatch", Bool_t tUseNewGeoData = false, Int_t tAnaType = doAnalysis::StdHZZllll);
		~HiggsAnalysis();

		// intializing the variables and tools
		int InitializeVar();
		void getPeriodEvent();	

		// cutFlow analysis
		int AnalyzeTree();
		int AnalyzeTreeEvent(int dataNumber);
		
		// outputs...
		void PrintInitVar();
		void printDebug(int EntNumber, Bool_t printllee = false);
		void printDebugInfo();
		void printMuonInfo();
		
		// Printing Functions
		void SetupPrintEventList(Bool_t overWrite, TString fileName);
		void PrintEventList(Bool_t passCut4Mu, Bool_t passCut4El, Bool_t passCut2L2L);
 		//void SetupPrintHist(Bool_t overWrite, TString fileName);		
 		void SaveHist(Bool_t overWrite);
		TString getSampleName();

		// For truth Studies
		void printMCInfo(int EntNumber);
		void printMCEventInfo(QuadLepton * higgs);	
		void printMCParticleInfo(Int_t index);
		
		TString getParticleName(int pdgID);
	private :
		// For Counting Events
		void FillCountingHist();

		// Functions for Channel Specific analysis
		Bool_t CutFlow4Mu(Double_t weight);
		Bool_t CutFlow4El(Double_t weight);
		Bool_t CutFlow2L2L(Double_t weight, Int_t forceQuadType = -1);
		Bool_t CutFlowllee(Double_t weight);

		// CutFlow Functions
		Bool_t DataPreselectionCut();
		Bool_t AllPreselectionCut();
		
		// Vertex Cut
		Bool_t VertexCut();
		Int_t getNVertex(Int_t nCutVertex);

		// BCH cut
		Bool_t BCHCut();

		// Trigger Functions
		Bool_t SingleElectronTrigger(Int_t runNumber);
		Bool_t DiElectronTrigger(Int_t runNumber);
		Bool_t SingleMuonTrigger(Int_t runNumber);
		Bool_t DiMuonTrigger();
		Bool_t ElectronMuonTrigger();
		
		// OverLap Function
		void RemoveOverlap();
		Bool_t ElectronMuonOverlap(D3PDReader::ElectronD3PDObjectElement *el_curr);
		Bool_t CaloMuonElectronOverlap(D3PDReader::MuonD3PDObjectElement *mu_curr);
		Bool_t JetsElectronOverlap(D3PDReader::JetD3PDObjectElement *jets_curr);

		// Combination
		// Dilepton
		vector<DiLepton *> GetDiLeptonComb(vector<ChargedLepton *> singleEvent);
		// Quad combinations
		vector<QuadLepton *> GetQuadLeptonComb(vector<DiLepton *> diLepton); // for 4Mu and 4E analysis
		vector<QuadLepton *> GetQuadLeptonComb(vector<DiLepton *> diMuLepton, vector<DiLepton *> diElLepton, Int_t forceQuadType = -1); // for 2Mu2El and 2El2Mu
		// Get a single Quad lepton
		QuadLepton* GetQuadEvent(vector<QuadLepton *> higgs, Int_t * cutEventPass, Int_t type, TH1D *cutPassHistW, Double_t weight);
		// Cuts on the single Quad lepton
		Bool_t CutQuadLepton(QuadLepton * higgs, Int_t * cutEventPass, Int_t type, TH1D *cutPassHistW, Double_t *weightInit, Int_t * prodCH, Bool_t dolleeFlow);	

		// Helper Functions for the Quad Cut
		Bool_t InterpolationZ2Cut(Double_t Z2Mass, Double_t QuadMass, Int_t type);
		Bool_t DeltaRCut(vector<ChargedLepton *> curr_lep);
		Bool_t JPsiVeto(QuadLepton * higgs, Int_t type);
		// Track Iso
		vector<Double_t> getTrackPT(vector<ChargedLepton *> curr_lep);
		Bool_t CutTrackIso(vector<ChargedLepton *> curr_lep, vector<Double_t> trackPT, Int_t type, QuadLepton* higgs);
		// Calo Iso
		vector<Double_t> getCaloET(vector<ChargedLepton *> curr_lep);	
		Double_t getCaloIsoCorrection(D3PDReader::ElectronD3PDObjectElement* el_curr);
		Bool_t CutCaloIso(vector<ChargedLepton *> curr_lep, vector<Double_t> caloPT, Int_t type, QuadLepton* higgs);
		// Impact Parameter significance
		Bool_t D0SigCut(vector<ChargedLepton *> curr_lep, Int_t type, QuadLepton* higgs);

		// Trigger Matching
		Bool_t TriggerMatch(vector<ChargedLepton *> curr_lep ,Int_t type);
		Bool_t TriggerMatchSingleMuon(vector<ChargedLepton *> curr_lep, TString triggerName[]);
		Bool_t TriggerMatchDiMuon(vector<ChargedLepton *> curr_lep, TString triggerName[]);
		Bool_t TriggerMatchSingleElectron(vector<ChargedLepton *> curr_lep, TString triggerName[]);
		Bool_t TriggerMatchDiElectron(vector<ChargedLepton *> curr_lep, TString triggerName[]);
		Bool_t TriggerMatchElectronMuon(vector<ChargedLepton *> curr_lep, TString triggerName[]);
		
		// Mass calculations 
		void MassCalc (QuadLepton * higgs, Int_t Type);		
		void CorrectFSR(QuadLepton * higgs, Int_t Type);
		void CorrectZMassConstraint(QuadLepton * higgs, Int_t muType = muonType::CB);
	
		// weights
		Double_t getEventWeight();
		Double_t getZVertexWeight();
		Double_t getPileupWeight();
		Double_t getggFWeight();
		Double_t getJHUWeight();				
		Double_t getTTBarVeto();
		Double_t getZZOverlapVeto();
		Double_t getZBBOverlapVeto();
		Double_t getHiggsWeight(QuadLepton * higgs);
		Double_t getTriggerWeight(QuadLepton * higgs);
		Double_t getLepEffWeight(QuadLepton * higgs);
		void fillHiggsWeight(QuadLepton * higgs);
		
		// CrossSection Weights
		Double_t getMCHiggsMass();
		void fillCrossSection(QuadLepton * higgs);
		Double_t getCrossSectionWeight();

		// Categorization
		void fillProductionChannel(QuadLepton * higgs, vector<ChargedLepton *> muObject, vector<ChargedLepton *> elObject, vector<ChargedLepton *> jetObjects);
		void fillBDTWeights(QuadLepton * higgs);
		Bool_t isGoodExtraLepton(QuadLepton * higgs, ChargedLepton * lep);
		void fillFudicialJets(QuadLepton * higgs);
		
		// For Stream decision and to take care of the overlap
		void fillStreamAnalysis();
		void fillFlagStream(QuadLepton * higgs);

		// Fill truth information
		void fillTruthRecoMatchedInfo(QuadLepton * higgs);
		Int_t getIndexBarcodeMatch(Int_t truthBarcode, Int_t truthPDG);
		Int_t getIndexDeltaRMatch(Int_t truthPDG, Double_t eta, Double_t phi);		
		void fillTruthInfo(QuadLepton * higgs, TLorentzVector bornSumVec);
		void fillTruthJetsInfo(QuadLepton * higgs);
		vector<Int_t> getBornIndexFromBare(vector<Int_t> lepBareTruthIndex, vector<Int_t> lepBareTruthPDGID, Bool_t isRecoMatched);
		vector<Int_t> getBornIndexFromBarePythia6(vector<Int_t> lepBareTruthIndex, vector<Int_t> lepBareTruthPDGID, Bool_t isRecoMatched);
		Bool_t checkSameParent(Int_t index1, Int_t index2);
		Bool_t containStatus3Taus();
		Int_t getTruthQuadType();
		Bool_t checkParentHiggs(Int_t index);
		Bool_t checkParent(Int_t index, Int_t parentPDGID, Int_t constParentIndex = -1);
		void fillExtraLepParent(ChargedLepton* lep);

		vector<Int_t> getMotherIndexFromBare(vector<Int_t> lepBareTruthIndex, vector<Int_t> lepBareTruthPDGID, vector<Int_t> lepBornTruthIndex);
		vector<TLorentzVector> getVecFromIndex(vector<Int_t> lepIndex);

		// LepId stuff
		void fillLepID(QuadLepton * higgs);
		void fillLepIDLepton(ChargedLepton * lep);

		// CR stuff
		vector<QuadLepton*> cutClosestZ1(vector<QuadLepton*> higgsContainer);
		void fillLooseCut(QuadLepton* higgs);

		// For getting muon ID and MS vars
		void fillMuonHelperVars();

		// Helper CutFlow Function
		void InitializeEventVar();
		void FillProductionTag();
		void InitPileupTool();
		void InitTriggerMatchingToolMain();
		void InitZVertexTool();		
		void InitTriggerMatchingTool();
		void InitGoodRunList();
		void InitTileTrip();
		void InitBCHCleaning();
		void InitggFReweight();
		void InitCategoryBDTTool();
		void InitBDTTool();
		void InitJHUReweight();
		void fillEventVarInfo(QuadLepton * higgs, Int_t type, Int_t *prodCH);
		
		
		void FillTriggerString(TString singleMu [], TString diMu[], TString singleEl [], TString diEl [], TString eMu[]);

		Double_t DeltaR (Double_t eta_1, Double_t phi_1, Double_t eta_2, Double_t phi_2);
		Double_t DeltaRCaloIso (ChargedLepton *lep1, ChargedLepton *lep2);
		
		void clearVar();

};

#endif


