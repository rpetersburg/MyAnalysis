#include <stdlib.h>
#include <string>
#include <math.h>
#include <iomanip>
#include <string>

#include "MyAnalysis/HiggsAnalysis.h"
#include "MyAnalysis/ChargedLepton.h"
#include "D3PDReader/MuonD3PDObject.h"
#include "D3PDReader/ElectronD3PDObject.h"
#include "MyAnalysis/macroDef.h"

#include <TTree.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <iostream>
#include <TDirectory.h>
#include <TString.h>
#include <TMath.h>
#include <TFile.h>
#include <TObjArray.h>


////////////////////////////////////////////////////////////////////////////////////////
//				Constructor and Destructor
////////////////////////////////////////////////////////////////////////////////////////
HiggsAnalysis::HiggsAnalysis (TTree *phyObject, Bool_t tRunningGrid, TString tgridFileName, Bool_t tUseNewGeoData, Int_t tAnaType)
{	

	// Reading the TTree
	event = new D3PDReader::Event();
	physicsTree = phyObject;
	event->ReadFrom(phyObject);

	//Initializing the Vars to zero for consistency
	curEvent = 0;
	dataYear = 0;
	CM_E = 0;
	curMCCollection = 0;
	dataPeriod = "";
	currFileName = "";
	sampleProdType = -1;
	noTauSample = false;
	anaType = tAnaType;
	generatorName = MCGeneratorName::other;
	countingHist = new TH1F("CountingHist", "CountingHist", 10, 0, 10);

	cout<<"--------------------------------"<<endl;
	if(anaType == doAnalysis::StdHZZllll) cout<<"H->ZZ->4l analysis"<<endl;
	else if(anaType == doAnalysis::Zllll) cout<<"Z->4l analysis"<<endl;
	else if(anaType == doAnalysis::trigeff4l) cout<<"H->ZZ->4l TrigEff Study"<<endl;
	cout<<"--------------------------------"<<endl<<endl;

	event->GetEntry(0);
	InitializeEventVar();
	InitializeVar();

	// For running on the grid
	// These vars will be specified by outside world if running on grid
	runningGrid = tRunningGrid;
	gridFileName = tgridFileName;
	cout<<"Constructor setup runningGrid: "<<runningGrid<<endl;
		
	// choosing the geo for the data
	useNewGeoData = tUseNewGeoData;
	cout<<"if Data, use newGeo?: "<<useNewGeoData<<endl;
	
	// Just a boolean that controls Trigger matchingh
	// Changes according to what data year it is
	D3DPTriggerDev = true;
	
	// This controls the debug printing. Gets set if the debug specific function is called
	isDebugCall = false;
	isDebugCallllee = false;
	printWeight = false;

	// What analysis to do CHOOSE HERE.
	do4Mu 	= true;
	do4El 	= true;
	do2L2L 	= true;
	dollee 	= false;
	doCR	= false;
	doSysTree = false;

	if(!doCR) dollee = false;
	// Overwrites the above decision based on stream name
	fillStreamAnalysis();

	// Setup the analysis
	doCorr 		= true;
	doWeight 	= true;

	doMoriondConfig = false;
	// For Setting up what electron ID tool to use. Either Likelihood(new) or Multilepton
	if(doMoriondConfig) useLikelihood = false;
	else useLikelihood = true;

	// So that I don't get to use any MC vars
	if(!isMC) doWeight = false;

	if(doWeight)
	{	
		// I can't have weights without correction... I think. I will have to check this
		doCorr = true;
		// For Pileup weight
		doPileupWeight = true;
		// For Zvertex
		//if(dataYear == 2012) doZVertexWeight = true;
		//else if(dataYear == 2011) doZVertexWeight = false;
		// New recommendation
		doZVertexWeight = true;
		// ggF Weight 
		if(dataYear == 2012) doggFWeight = false;
		else if(dataYear == 2011) doggFWeight = true;
		// JHU weight
		doJHUWeight = true;

		// For TTbar overlap
		doTTBarOverlap = true; // making it false for now.. makes the program too slow.. I Think
		// For ZZOverlap
		doZZOverlap = true; // making it false for now.. makes the program too slow .. I think
		// For ZBBOverlap
		doZBBOverlap = true; // making it false for now.. makes the program too slow .. I think
		// ID and Reco Efficiency
		doScaleEfficiency  = true;
		// Trigger Effciency SF
		doTriggerEfficiency = true;
	}
	else
	{
		doPileupWeight = false;
		doZVertexWeight = false;
		doggFWeight = false;
		doJHUWeight = false;
		doTTBarOverlap = false;
		doZZOverlap = false;
		doZBBOverlap = false;
	}
	// Setting up the corrections
	if(doCorr)
	{	
		// D0/Z0 smearing setup based on year
		if(dataYear == 2012) doSmearD0Z0 = true;
		else if(dataYear == 2011) doSmearD0Z0 = false;
		else cout<<"Smear D0Z0: Year not recognized"<<endl;
		
		// MC Smear
		doSmearMC = true;

		// Data Scaling
		doScaleData = true;
		// Muon Scaling for uncertainities
		doScaleMuon = true;
		
		// scale crak setup based on year
		if(dataYear == 2012) doScaleCrack = false;
		else if(dataYear == 2011) doScaleCrack = true;
		else cout<<"Scale Crack: Year not recognized"<<endl;
		
		// Jet calibration smearsetup
		doJetCal = true;	
		
		// Calo Iso Correction
		doCaloIsolationCorr = true;

		// Electron cluster and track combination
		doEPcombination = true;

	}
	else
	{
		doSmearD0Z0 = false;
		doSmearMC = false;
		doScaleMuon = false;
		doScaleData = false;
		doScaleCrack = false;
		doCaloIsolationCorr = false;
		doJetCal = false;
		doEPcombination = false;
	}
	
	// Setting the lepton Collection
	if(dataYear == 2011) electronCollection = electronCollection::LoosePlusPlus;
	else if(dataYear == 2012 && !useLikelihood) electronCollection = electronCollection::MultiLepton;
	else if(dataYear == 2012 && useLikelihood) electronCollection = electronCollection::Likelihood;
	else cout<<"Error LeptonCollection: dataYear not recongnized"<<endl;
	muonCollection = muonCollection::Loose;
	// For corrected calo iso
	if(dataYear == 2012)
	{
		CaloIsoCorrection::SetPtLeakageCorrectionsFile("../../egammaAnalysisUtils/share/isolation_leakage_corrections.root");
	}

	TChain* chain = dynamic_cast<TChain *> (physicsTree);	
	currFileName = chain->GetFile()->GetPath();	
	event->GetEntry(0);
	InitializeEventVar();
	InitializeVar();


	// calling the initializing functions to setup the above vars
	InitPileupTool();
	InitZVertexTool();
	InitGoodRunList();
	InitTileTrip();
	InitCategoryBDTTool();
	InitBDTTool();
	InitJHUReweight();
	InitggFReweight();	
	higgs_bror = new H4lBrRatio();

	// Correction
	// Init the class that contains all the function
	corr = new Correction(dataYear, isMC, doSmearMC, doScaleData, doScaleCrack, 
			doJetCal, doScaleMuon, doTriggerEfficiency,doMoriondConfig,
			curMCCollection, curDataCollection, doEPcombination, pileupTool);
	corr->InitD0Z0Smear();
	corr->InitMuonSmear();
	corr->InitElectronSmear(electronCollection);
	corr->InitJetCal();

	// Initializing the Hist
	Hist = new HistContainer(nCut, nMuCut, nElCut, nJetsCut, nCH);

	// Init the lepton tools
	electronObj = new ElectronObject(dataYear, useLikelihood);
	electronObjLoose = new ElectronObject(dataYear, useLikelihood);
	
	Bool_t doZ4lLooseCut = false;
	if(anaType == doAnalysis::Zllll) doZ4lLooseCut = true;
	muonObj = new MuonObject(dataYear, doZ4lLooseCut);
	
	jetsObj = new JetsObject(dataYear);
	jetsObjTruth = new JetsObject(dataYear);
	jetsObj_Fid = new JetsObject(dataYear);
	jetsObjTruth_Fid = new JetsObject(dataYear);

	// Output Tree
	outputTree = new OutputTree();
	outputFilePath = "Output/output_giveFileNamePlease.root"; // Just in case the outside world didn't specify it
	
	outputTreeCR = new OutputTree();
	outputTreelleeCR = new OutputTree();

	if(doSysTree) outputTreeSys = new OutputTreeSys();
	
	
	// Initializing trigger matching
    InitTriggerMatchingToolMain(); 

	// Print event list
	printEventList = true;
	printMass = true;

	// Setting up the type of calibration that we are doing
	if (!doCorr) curCalibration = calibrationType::noCalib;
	else{
		if (curMCCollection == MCCollection::MC11c ||
			curMCCollection == MCCollection::MC12a ||
			curMCCollection == MCCollection::MC12b ||
			curDataCollection == dataCalibType::y2011c ||
			curDataCollection == dataCalibType::y2012ab ){
			if (doEPcombination)	curCalibration = calibrationType::stdCalibEp;
			else 					curCalibration = calibrationType::stdCalib;
		}
		else if (curMCCollection == MCCollection::MC11d ||
				 curMCCollection == MCCollection::MC12c ||
			curDataCollection == dataCalibType::y2011d ||
			curDataCollection == dataCalibType::y2012c ){
			if (doEPcombination)	curCalibration=calibrationType::MvaCalibEp;
			else					curCalibration=calibrationType::MvaCalib;
		}
	}
	
}

HiggsAnalysis::~HiggsAnalysis ()
{
	//return;
	// Don't delete if it is a debug call. Something goes wrong if you do
	//if(isDebugCall) return;
	delete event;
	delete muonObj;
	delete electronObj;
	delete electronObjLoose;
	delete jetsObj;
	delete jetsObjTruth;
	delete jetsObj_Fid;
	delete jetsObjTruth_Fid;
	delete outputTree;
	delete outputTreeCR;
	delete outputTreelleeCR;
	if(doSysTree) delete outputTreeSys;
	//delete higgsCandidate4Mu;
	//delete higgsCandidate4El;
	//delete higgsCandidate2L2L;

	delete[] cutPass;
	delete[] cutMuPass;
	delete[] cutElPass;
	delete[] cutJetsPass;
	delete[] cutPassW;
	delete[] cut4MuPass;
	delete[] cut4ElPass;
	delete[] cut2L2LPass;
	delete[] cut4MuPassW;
	delete[] cut4ElPassW;
	delete[] cut2L2LPassW;
	delete[] prodCH4Mu;
	delete[] prodCH4El;
	delete[] prodCH2L2L;
	delete[] truthQuadType;
	delete[] nTruthQuadType;
	delete countingHist;


	// Tools
	delete pileupTool;
	delete triggerNavigationVariables;
	delete muonTriggerMatchTool;
	delete electronTriggerMatchTool;
	delete leptonSF;
	delete zVertexTool;
	delete TileTrip;
	delete higgs_bror;
	delete ggFReweight;
	delete corr;
	delete CategoriesDiscriminantTool;   
	delete Hist;
	delete ptJHUReweighting;
	delete BDTtool;
}


////////////////////////////////////////////////////////////////////////////////////////
//				Analysis Functions
////////////////////////////////////////////////////////////////////////////////////////
// Analyzes the whole TTree by calling the analyzeTreeEvent on each event
int HiggsAnalysis::AnalyzeTree()
{
	
	// For looping over the trees
	Int_t currEvent = 0;
	// To count the events
	Int_t countPassed = 0;
	
	// For the main Loop
	Int_t nEvent = physicsTree->GetEntries();
	// Main loop
	for(Long64_t iEvent = 0; iEvent < nEvent; iEvent++)
	{
		curEvent = iEvent;
		Long64_t currEvent = iEvent;
		TChain* chain = dynamic_cast<TChain *> (physicsTree);
		if(chain)
		{
			currEvent = chain->LoadTree(currEvent);
		}
		currFileName = chain->GetFile()->GetPath();
		Double_t mHiggs = getMCHiggsMass();
		if(AnalyzeTreeEvent(currEvent)) countPassed++;
		
		// For keeping track
		if(iEvent % 5000 == 0) cout<<"Current Event: "<<iEvent<<endl;
	}
	
	// Fill the counting Hist
	FillCountingHist();
	// Closing the file
	file4Mu.close();
	file4El.close();
	file2El2Mu.close();
	file2Mu2El.close();
	file2L2L.close();
	fileSampleName.close();

	// Saving the tree
	outputTree->saveTrees(outputFilePath, countingHist, getSampleName());
	
	TString outputFilePathCR = outputFilePath(0, outputFilePath.Length()-5) + "_CR.root";
	outputTreeCR->saveTrees(outputFilePathCR, countingHist, getSampleName());
	
	TString outputFilePathlleeCR = outputFilePath(0, outputFilePath.Length()-5) + "_lleeCR.root";
	outputTreelleeCR->saveTrees(outputFilePathlleeCR, countingHist, getSampleName());

	if(doSysTree)
	{
		TString outputFilePathSys = outputFilePath(0, outputFilePath.Length()-5) + "_Sys.root";
		outputTreeSys->saveTrees(outputFilePathSys, countingHist, getSampleName());
	}
	return countPassed;
}

// Analyze each event
int HiggsAnalysis::AnalyzeTreeEvent(int dataNumber)
{
	// Loading the event
	if(dataNumber >= 0) event->GetEntry(dataNumber);

	// Getting the truth quad type
	nTruthQuadType[getTruthQuadType()]++;
	Hist->truthQuadHist->Fill(getTruthQuadType(), 1);

	// Printing the currentMC collection
	if(isDebugCall || curEvent == 0)
	{
		cout << fixed;		
		cout<<"--------------------------------------"<<endl;
		if(curMCCollection == MCCollection::MC11c) cout<<"MC11c collection"<<endl;
		else if(curMCCollection == MCCollection::MC12a) cout<<"MC12a collection"<<endl;
		else if(curMCCollection == MCCollection::MC12b) cout<<"MC12b collection"<<endl;
		else if(curMCCollection == MCCollection::MC12c) cout<<"MC12c collection"<<endl;
		else if(curDataCollection == dataCalibType::y2011c) cout<<"data_11 Old Geo"<<endl;
		else if(curDataCollection == dataCalibType::y2011d) cout<<"data_11 New Geo"<<endl;
		else if(curDataCollection == dataCalibType::y2012ab) cout<<"data_12 Old Geo"<<endl;
		else if(curDataCollection == dataCalibType::y2012c) cout<<"data_12 New Geo"<<endl;
		cout<<"--------------------------------------"<<endl;		
	}

	// Getting the initial event weight
	Double_t eventWeight = 1.;
	if(doWeight && isMC) eventWeight = eventWeight * getEventWeight();
	
	// Filling the histogram for counting
	countingHist->Fill(1);
	Double_t eventWeight_raw = eventWeight;
	eventWeight_raw = eventWeight_raw/getJHUWeight();
	countingHist->Fill(2,eventWeight); // To get rid of the ggF Weight
	countingHist->Fill(3,eventWeight);
	countingHist->Fill(4,(eventWeight/getggFWeight())/getJHUWeight());
	countingHist->Fill(5,eventWeight/getJHUWeight());
	// Setting the weight for the plot histograms
	Hist->weight = eventWeight * getCrossSectionWeight();

	cutPass[cutFlow::Total] ++;
	cutMuPass[cutMuFlow::Total] += (event->mu_staco.n() + event->mu_calo.n());
	cutElPass[cutElFlow::Total] += (el_cur->n());
	cutJetsPass[cutJetsFlow::Total] += event->jet_akt4topoem.n();
	Hist->cutPassHistW->Fill(cutFlow::Total, eventWeight);

	// Initializing the event Specific Variables
	InitializeEventVar();
	getPeriodEvent();
	Bool_t passCut = true;
	
	// Preselection cut for Data
	if(!isMC)
	{
		passCut = DataPreselectionCut();
		if(!passCut) return passCut;
	}
	cutPass[cutFlow::DataPreselection]++;
	cutMuPass[cutMuFlow::DataPreselection] += (event->mu_staco.n() + event->mu_calo.n());
	cutElPass[cutElFlow::DataPreselection] += (el_cur->n());
	cutJetsPass[cutJetsFlow::DataPreselection] += event->jet_akt4topoem.n();
	Hist->cutPassHistW->Fill(cutFlow::DataPreselection, eventWeight);	
	
	//Rest of Preselection
	passCut = AllPreselectionCut();
	if(!passCut) return passCut;
	cutPass[cutFlow::Preselection]++;
	cutMuPass[cutMuFlow::Preselection] += (event->mu_staco.n() + event->mu_calo.n());
	cutElPass[cutElFlow::Preselection] += (el_cur->n());
	cutJetsPass[cutJetsFlow::Preselection] += event->jet_akt4topoem.n();
	Hist->cutPassHistW->Fill(cutFlow::Preselection, eventWeight);

	// Initial Trigger Cut
	if(isMC)
	{ 
		pileupTool->SetRandomSeed(314159+event->eventinfo.mc_channel_number()*2718+event->eventinfo.EventNumber());
		runNumber_sf = pileupTool->GetRandomRunNumber(event->eventinfo.RunNumber());
	}
	else runNumber_sf = event->eventinfo.RunNumber();
	
	// For Counting DiMuon Trigger
	Bool_t passCut4Mu = SingleMuonTrigger(runNumber_sf) | DiMuonTrigger() ;
	if(passCut4Mu) {cutPass[cutFlow::Trigger4Mu] ++; Hist->cutPassHistW->Fill(cutFlow::Trigger4Mu, eventWeight);}
	
	// for Counting electron trigger
	Bool_t passCut4e = SingleElectronTrigger(runNumber_sf)| DiElectronTrigger(runNumber_sf);
	if(passCut4e) {cutPass[cutFlow::Trigger4e] ++;; Hist->cutPassHistW->Fill(cutFlow::Trigger4e, eventWeight);}
	
	// If it passes any of the trigger, move ahead
	passCut = SingleElectronTrigger(runNumber_sf)| DiElectronTrigger(runNumber_sf) | 
	SingleMuonTrigger(runNumber_sf) | DiMuonTrigger() |
	ElectronMuonTrigger();
	
	// For Trig eff studies
	if(anaType == doAnalysis::trigeff4l) passCut = true;

	if(!passCut) return passCut;
	
	cutPass[cutFlow::Trigger]++;
	cutMuPass[cutMuFlow::Trigger] += (event->mu_staco.n() + event->mu_calo.n());
	cutElPass[cutElFlow::Trigger] += (el_cur->n());
	cutJetsPass[cutJetsFlow::Trigger] += event->jet_akt4topoem.n();
	Hist->cutPassHistW->Fill(cutFlow::Trigger, eventWeight);

	// Vars to store the efficiency
	vector<Double_t> muonStacoEff;
	vector<Double_t> muonCaloEff;

	// Performs the Cuts on Muon..
	// Clearing the vars
	muonObj->clearVars();	
	muonObj->SetHist(Hist);

	// Clearing the corr
	corr->ClearVars();
	
	fillMuonHelperVars();

	// Smearing
	if(isDebugCall) printMuonInfo();
	if(isDebugCall) corr->debugCall();
	if(doCorr && doSmearD0Z0 && isMC) corr->SmearD0Z0(&(event->mu_staco), event->eventinfo.EventNumber(), leptonType::MuonStaco);
	if(doCorr && doSmearD0Z0 && isMC) corr->SmearD0Z0(&(event->mu_calo), event->eventinfo.EventNumber(), leptonType::MuonCalo);
	if((doCorr && doSmearMC) || doScaleData) corr->SmearMuon(&(event->mu_calo), event->eventinfo.EventNumber(), leptonType::MuonCalo);
	if((doCorr && doSmearMC) || doScaleData) corr->SmearMuon(&(event->mu_staco), event->eventinfo.EventNumber(), leptonType::MuonStaco);
	if(isDebugCall) printMuonInfo();

	// getting the Eff values
	if(doWeight && doScaleEfficiency) 
	{
		muonStacoEff = corr->muonStacoEff;
		muonCaloEff = corr->muonCaloEff;
	}
	// Otherwise just fill them with 1
	else
	{
		for(Int_t i = 0; i < event->mu_staco.n(); i++)
			muonStacoEff.push_back(1);

		for(Int_t i = 0; i < event->mu_calo.n(); i++)
			muonCaloEff.push_back(1);
	}
	
	muonObj->FillMuon(&(event->mu_staco), leptonType::MuonStaco, muonStacoEff, isMC);
	muonObj->FillMuon(&(event->mu_calo), leptonType::MuonCalo, muonCaloEff, isMC);

	Bool_t passCutMu = muonObj->MuonCut(cutMuPass);
	if(!passCut) return passCut;

	// To Perform cut on electrons
	// Var to store the smearing
	vector<Double_t> elSmearVal;
	vector<Double_t> elEff;
	vector<Double_t> elResolution;

	// Clearing the vars
	electronObj->clearVars();	
	electronObj->SetHist(Hist);	
	electronObjLoose->clearVars();	
	electronObjLoose->SetHist(Hist);	
	// Smearing
	if(doCorr && doSmearD0Z0 && isMC) corr->SmearD0Z0(el_cur, event->eventinfo.EventNumber());
	if((doCorr && doSmearMC) || doScaleData || doEPcombination) 
	{
		
		// Don't want to correct Nom in this loop. This will create a double smearing
		if(doSysTree)
		{
			for(Int_t i = 0; i < doSys::Nom; i++)
			{
				elSmearVal = corr->SmearElectron(el_cur, event->eventinfo.EventNumber(), runNumber_sf, i);
				corr->ClearVars();
			}
		}
		corr->ClearVars();	
		elSmearVal.clear();		
		elSmearVal = corr->SmearElectron(el_cur, event->eventinfo.EventNumber(), runNumber_sf);
	}
	// In cases where no smearing will be done
	if(!doCorr || !doSmearMC) 
	{
		elSmearVal.clear();
		for(Int_t i = 0; i < el_cur->n(); i++)
			elSmearVal.push_back(1);
	}
	// getting the Eff values
	if(doWeight && doScaleEfficiency) elEff = corr->electronEff;
	// Otherwise just fill them with 1
	else 
	{
		elEff.clear();
		for(Int_t i = 0; i < el_cur->n(); i++)
			elEff.push_back(1);
	}
	elResolution = corr->electronEpErr;
	if(elResolution.size() == 0)
	{
		for(Int_t i = 0; i < el_cur->n(); i++)
			elResolution.push_back(-1);
	}
	// Getting the clusterPt for getting the lepton ID
	vector<Double_t> el_bfEP_clPt;
	el_bfEP_clPt.clear();
	// Sanity check
 	if( el_cur->n() != (Int_t) corr->bfEP_cl_pt.size() && doCorr) cout<<"Error: beforeEP_cl_pt has a different size than electron container"<<endl;
 	for(Int_t i = 0; i < el_cur->n(); i++)
 	{
 		if(doCorr) el_bfEP_clPt.push_back(corr->bfEP_cl_pt[i]);
 		else el_bfEP_clPt.push_back((*el_cur)[i].cl_pt());
 	}

	electronObj->FillElectron(el_cur, leptonType::ElectronGSF, elSmearVal, elEff,elResolution, el_bfEP_clPt ,isMC);
	electronObjLoose->FillElectron(el_cur, leptonType::ElectronGSF, elSmearVal, elEff,elResolution, el_bfEP_clPt, isMC);

	Bool_t passCutEl = electronObj->ElectronCut(cutElPass, getNVertex(2));
	passCutEl = electronObjLoose->ElectronCut(cutElLoosePass, getNVertex(2), true);

	// Photon Smearing
	vector<Double_t> phSmearVal;
	if((doCorr && doSmearMC) || doScaleData) phSmearVal = corr->SmearPhoton(&(event->ph), event->eventinfo.EventNumber(), runNumber_sf);
	// In cases where no smearing will be done
	if(!doCorr || !doSmearMC) 
	{
		phSmearVal.clear();
		for(Int_t i = 0; i < el_cur->n(); i++)
			phSmearVal.push_back(1);
	}

	// To Perform cut on jets 
	// Clearing the vars
	jetsObj->clearVars();
	jetsObjTruth->clearVars();
	jetsObj_Fid->clearVars();
	jetsObjTruth_Fid->clearVars();
	// Calibrate jets
	Double_t mu = -9999;
	mu = event->eventinfo.averageIntPerXing();
	if(isMC)
	{
		if(dataYear == 2012)  mu = (event->eventinfo.lbn()==1&&int(mu+0.5)==1)?0.:mu;
	}
	Hist->muInteractionHist->Fill(mu);
	Hist->muInteractionHistPileupW->Fill(mu,getPileupWeight());

	Double_t rhoKt = -1;
	if(dataYear == 2012) rhoKt = event->Eventshape.rhoKt4EM();
	if(doCorr && doJetCal) corr->CalibrateJet(&(event->jet_akt4topoem), dataYear, rhoKt, mu, getNVertex(2));
	jetsObj->FillJets(&(event->jet_akt4topoem), jetsType::AntiKt4TopoEM, isMC, event->eventinfo.EventNumber());
	Bool_t passCutJets = jetsObj->JetsCut(cutJetsPass, event->eventinfo.RunNumber());
	// True jets
	jetsObjTruth->FillJetsTruth(&(event->jet_antikt4truth), jetsType::AntiKt4TopoEMTruth, isMC);
	Bool_t passCutJetsTruth = jetsObjTruth->JetsCut(cutJetsPass, event->eventinfo.RunNumber());

	// Jets for Fudical Cross-Section
	jetsObj_Fid->FillJets(&(event->jet_akt4topoem), jetsType::AntiKt4TopoEM_Fid, isMC, event->eventinfo.EventNumber());
	Bool_t passCutJets_Fid = jetsObj_Fid->JetsCut(cutJetsPass, event->eventinfo.RunNumber());
	// True jets
	jetsObjTruth_Fid->FillJetsTruth(&(event->jet_antikt4truth), jetsType::AntiKt4TopoEMTruth_Fid, isMC);
	Bool_t passCutJetsTruth_Fid = jetsObjTruth_Fid->JetsCut(cutJetsPass, event->eventinfo.RunNumber());

	// Reading the muon, electrons and jet events
	clearVar();
	jetsOverlap = jetsObj->getJetsVec();	
	muOverlap = muonObj->getMuonVec();
	elOverlap = electronObj->getElectronVec();
	elLooseOverlap = electronObjLoose->getElectronVec();
	jetsTruthEvent = jetsObjTruth->getJetsVec();
	jetsOverlap_Fid = jetsObj_Fid->getJetsVec();
	jetsTruthEvent_Fid = jetsObjTruth_Fid->getJetsVec();

	// Removing the overlap
	RemoveOverlap();

	// Channel Specific Cuts
	Bool_t pass4MuCut = false;
	Bool_t pass4ElCut = false;
	Bool_t pass2L2LCut = false;
	Bool_t passlleeCut = false;
	Bool_t passCutCH = false;
	if(do4Mu && !passCutCH) pass4MuCut = CutFlow4Mu(eventWeight);
	passCutCH = passCutCH | pass4MuCut;

	if(do2L2L && !passCutCH) pass2L2LCut = CutFlow2L2L(eventWeight, quadType::El2Mu2);
	passCutCH = passCutCH | pass2L2LCut;

	if(do2L2L && !passCutCH) pass2L2LCut = CutFlow2L2L(eventWeight, quadType::Mu2El2);
	passCutCH = passCutCH | pass2L2LCut;

	//if(do2L2L && !passCutCH) pass2L2LCut = CutFlow2L2L(eventWeight);
	//passCutCH = pass2L2LCut;

	if(do4El && !passCutCH) pass4ElCut = CutFlow4El(eventWeight);
	passCutCH = passCutCH | pass4ElCut;	

	passCut = pass4MuCut | pass4ElCut | pass2L2LCut;
	if(pass4MuCut) {cut4MuPass[cutFlowCH::Final]++; Hist->cut4MuPassHistW->Fill(cutFlowCH::Final);}
	if(pass4ElCut) {cut4ElPass[cutFlowCH::Final]++; Hist->cut4ElPassHistW->Fill(cutFlowCH::Final);}
	if(pass2L2LCut) {cut2L2LPass[cutFlowCH::Final]++; Hist->cut2L2LPassHistW->Fill(cutFlowCH::Final);}


	// CR
	if(dollee) passlleeCut = CutFlowllee(eventWeight);

	// Printing the event list
	if(printEventList) PrintEventList(pass4MuCut, pass4ElCut, pass2L2LCut);
	
	// Cleaning up for this event
	muonStacoEff.clear();	
	muonCaloEff.clear();
	elSmearVal.clear();
	elEff.clear();
	elResolution.clear();
	phSmearVal.clear();


	return passCut;
}
////////////////////////////////////////////////////////////////////////////////////////
// 			Channel specific cutflow
//		all Function return 1 or true if event passed the cut
////////////////////////////////////////////////////////////////////////////////////////
Bool_t HiggsAnalysis::CutFlow4Mu(Double_t weight)
{
	// Counting the initial events
	cut4MuPass[cutFlowCH::Total]++;
	Hist->cut4MuPassHistW->Fill(cutFlowCH::Total, weight);

	// Imposing Channel specific Trigger
	Bool_t passCut4Mu = SingleMuonTrigger(runNumber_sf) | DiMuonTrigger();
	// For Trig eff studies
	if(anaType == doAnalysis::trigeff4l) passCut4Mu = true;
	
	if(!passCut4Mu) return false;
	cut4MuPass[cutFlowCH::Trigger]++;
	Hist->cut4MuPassHistW->Fill(cutFlowCH::Trigger, weight);

	// Ensure that event has the minimum number of required leptons: 4mu
	if(!(muEvent.size() >= 4)) return false;
	cut4MuPass[cutFlowCH::Lepton]++;
	Hist->cut4MuPassHistW->Fill(cutFlowCH::Lepton, weight);

	// Getting the diLepton combination
	diEvent2Mu.clear();
	diEvent2Mu = GetDiLeptonComb(muEvent);
	eventToCleanDi.insert(eventToCleanDi.end(), diEvent2Mu.begin(), diEvent2Mu.end());
	if(!(diEvent2Mu.size() >= 2)) return false;
	
	// Get the quadletpon combinations
	vector<QuadLepton *> higgsContainer =  GetQuadLeptonComb(diEvent2Mu);
	eventToClean4mu = higgsContainer;
	if(higgsContainer.size()<1) return false;
	cut4MuPass[cutFlowCH::SFOS]++;
	Hist->cut4MuPassHistW->Fill(cutFlowCH::SFOS, weight);
	
	// Get a single Quad Event based on SFOS and PGD mass
	QuadLepton* higgs = GetQuadEvent(higgsContainer, cut4MuPass, analysisType::Mu4, Hist->cut4MuPassHistW, weight);
	if(higgs == 0) return false;

	// Selection cuts on the the quadrilepton
	Bool_t passCut = CutQuadLepton(higgs, cut4MuPass, analysisType::Mu4, Hist->cut4MuPassHistW, &weight, prodCH4Mu, false);
	if(!passCut) return passCut;
	
	// Choosing the Final Higgs Candidate
	higgsCandidate4Mu = higgs;
	
	// Saving the weight
	higgsCandidate4Mu->weight_corr = weight;
	
	//fillEventVarInfo(higgsCandidate4Mu, analysisType::Mu4, prodCH4Mu);

	// Filling Histogram
	Hist->hist4MuMUnconstrained->Fill(higgsCandidate4Mu->getMass(), Hist->weight);
	Hist->hist4MuMFSR->Fill(higgsCandidate4Mu->getMassFSR(), Hist->weight);
	Hist->hist4MuMConstrained->Fill(higgsCandidate4Mu->getMassZMassCons(), Hist->weight);

	// Filling the tree
	//outputTree->fillTree(event, higgsCandidate4Mu, analysisType::Mu4, isMC);

	//printMCEventInfo(higgsCandidate4Mu);
	return true;
}

Bool_t HiggsAnalysis::CutFlow4El(Double_t weight)
{

	// Counting the inital events
	cut4ElPass[cutFlowCH::Total]++;
	Hist->cut4ElPassHistW->Fill(cutFlowCH::Total, weight);
	
	// Imposing Channel specific Trigger
	Bool_t passCut4e = SingleElectronTrigger(runNumber_sf)| DiElectronTrigger(runNumber_sf);
	// For Trig eff studies
	if(anaType == doAnalysis::trigeff4l) passCut4e = true;

	if(!passCut4e) return false;
	cut4ElPass[cutFlowCH::Trigger]++;
	Hist->cut4ElPassHistW->Fill(cutFlowCH::Trigger, weight);

	// Ensure that event has the minimum number of required leptons: 4e
	if(!(elEvent.size() >= 4)) return false;
	cut4ElPass[cutFlowCH::Lepton]++;
	Hist->cut4ElPassHistW->Fill(cutFlowCH::Lepton, weight);

	// Getting the diLepton combination
	diEvent2El.clear();
	diEvent2El = GetDiLeptonComb(elEvent);
	eventToCleanDi.insert(eventToCleanDi.end(), diEvent2El.begin(), diEvent2El.end());	
	if(!(diEvent2El.size() >= 2)) return false;

	// Get the quadletpon combinations
	vector<QuadLepton *> higgsContainer =  GetQuadLeptonComb(diEvent2El);
	eventToClean4el = higgsContainer;	
	if(higgsContainer.size()<1) return false;
	cut4ElPass[cutFlowCH::SFOS]++;
	Hist->cut4ElPassHistW->Fill(cutFlowCH::SFOS, weight);

	//cout<<"EventNumber: "<<event->eventinfo.EventNumber()<<endl;
	// Get a single Quad Event based on SFOS and PGD mass
	QuadLepton* higgs = GetQuadEvent(higgsContainer, cut4ElPass, analysisType::El4, Hist->cut4ElPassHistW, weight);
	if(higgs == 0) return false;

	// Selection cuts on the the quadrilepton
	Bool_t passCut = CutQuadLepton(higgs, cut4ElPass, analysisType::El4, Hist->cut4ElPassHistW, &weight, prodCH4El, false);
	if(!passCut) return passCut;
	
	// Choosing the Final Higgs Candidate
	higgsCandidate4El = higgs;

	// Saving the weight
	higgsCandidate4El->weight_corr = weight;
	
	//fillEventVarInfo(higgsCandidate4El, analysisType::El4, prodCH4El);
	
	// Filling Histrogram
	Hist->hist4ElMUnconstrained->Fill(higgsCandidate4El->getMass(), Hist->weight);
	Hist->hist4ElMFSR->Fill(higgsCandidate4El->getMassFSR(), Hist->weight);
	Hist->hist4ElMConstrained->Fill(higgsCandidate4El->getMassZMassCons(), Hist->weight);
	
	// Filling the tree
	//outputTree->fillTree(event, higgsCandidate4El, analysisType::El4, isMC);
	
	//printMCEventInfo(higgsCandidate4El);
	
	return true;
}


Bool_t HiggsAnalysis::CutFlow2L2L(Double_t weight, Int_t forceQuadType)
{
	// Counting the inital events
	cut2L2LPass[cutFlowCH::Total]++;
	Hist->cut2L2LPassHistW->Fill(cutFlowCH::Total, weight);

	// Imposing Channel specific Trigger
	Bool_t passCut2El2Mu = SingleElectronTrigger(runNumber_sf)| DiElectronTrigger(runNumber_sf) | 
	SingleMuonTrigger(runNumber_sf) | DiMuonTrigger() |
	ElectronMuonTrigger();
	// For Trig eff studies
	if(anaType == doAnalysis::trigeff4l) passCut2El2Mu = true;

	if(!passCut2El2Mu) return false;
	cut2L2LPass[cutFlowCH::Trigger]++;
	Hist->cut2L2LPassHistW->Fill(cutFlowCH::Trigger, weight);

	// Ensure that event has the minimum number of required leptons: 2e and 2mu
	if(!(elEvent.size() >= 2 &&  muEvent.size() >= 2 )) return false;
	cut2L2LPass[cutFlowCH::Lepton]++;
	Hist->cut2L2LPassHistW->Fill(cutFlowCH::Lepton, weight);
	
	// Getting the diLepton combination
	diEvent2El.clear();
	diEvent2El = GetDiLeptonComb(elEvent);
	diEvent2Mu.clear();
	diEvent2Mu = GetDiLeptonComb(muEvent);
	eventToCleanDi.insert(eventToCleanDi.end(), diEvent2El.begin(), diEvent2El.end());	
	eventToCleanDi.insert(eventToCleanDi.end(), diEvent2Mu.begin(), diEvent2Mu.end());	
	
	if(!(diEvent2El.size() >= 1 && diEvent2Mu.size() >= 1)) return false;
	
	// Get a single Quad Event based on SFOS and PGD mass
	vector<QuadLepton *> higgsContainer =  GetQuadLeptonComb(diEvent2Mu, diEvent2El, forceQuadType);
	eventToClean2mu2el = higgsContainer;	
	if(higgsContainer.size()<1) return false;
	cut2L2LPass[cutFlowCH::SFOS]++;
	Hist->cut2L2LPassHistW->Fill(cutFlowCH::SFOS, weight);

	// Getting the quadlepton
	QuadLepton* higgs = GetQuadEvent(higgsContainer, cut2L2LPass, analysisType::El2Mu2, Hist->cut2L2LPassHistW, weight);
	if(higgs == 0) return false;

	// Change the type depending on the current Higgs Candidate;
	Int_t type = analysisType::Mu2El2;
	if(higgs->getQuadType() == quadType::Mu2El2) {type = analysisType::Mu2El2;}
	else if(higgs->getQuadType() == quadType::El2Mu2) {type = analysisType::El2Mu2;}
	else {cout<<"ERROR: CutFlow2L2L: Type not recognized"<<endl;}

	// Selection cuts on the the quadrilepton
	Bool_t passCut = CutQuadLepton(higgs, cut2L2LPass, type, Hist->cut2L2LPassHistW, &weight, prodCH2L2L, false);
	if(!passCut) return passCut;
	// Choosing the Final Higgs Candidate
	higgsCandidate2L2L = higgs;

	// Saving the weight
	higgsCandidate2L2L->weight_corr = weight;
	
	//fillEventVarInfo(higgsCandidate2L2L, type, prodCH2L2L);
	
	// Filling Histrogram	
	if(type == analysisType::El2Mu2)
	{
		Hist->hist2El2MuMUnconstrained->Fill(higgsCandidate2L2L->getMass(), Hist->weight);
		Hist->hist2El2MuMFSR->Fill(higgsCandidate2L2L->getMassFSR(), Hist->weight);
		Hist->hist2El2MuMConstrained->Fill(higgsCandidate2L2L->getMassZMassCons(), Hist->weight);
	}
	else if(type == analysisType::Mu2El2)
	{
		Hist->hist2Mu2ElMUnconstrained->Fill(higgsCandidate2L2L->getMass(), Hist->weight);
		Hist->hist2Mu2ElMFSR->Fill(higgsCandidate2L2L->getMassFSR(), Hist->weight);
		Hist->hist2Mu2ElMConstrained->Fill(higgsCandidate2L2L->getMassZMassCons(), Hist->weight);
	}

	// Filling the tree
	//outputTree->fillTree(event, higgsCandidate2L2L, type, isMC);
	
	//printMCEventInfo(higgsCandidate2L2L);
	
	return true;
}

Bool_t HiggsAnalysis::CutFlowllee(Double_t weight)
{
	if(isDebugCall) isDebugCall = false; 
	if(isDebugCallllee) isDebugCall = true;
	// Counting the inital events
	cutlleePass[cutFlowCH::Total]++;
	Hist->cutlleePassHistW->Fill(cutFlowCH::Total, weight);

	// Imposing Channel specific Trigger
	Bool_t passCut2El2Mu = SingleElectronTrigger(runNumber_sf)| DiElectronTrigger(runNumber_sf) | 
	SingleMuonTrigger(runNumber_sf) | DiMuonTrigger() |
	ElectronMuonTrigger();
	// For Trig eff studies
	if(anaType == doAnalysis::trigeff4l) passCut2El2Mu = true;

	if(!passCut2El2Mu) return false;
	cutlleePass[cutFlowCH::Trigger]++;
	Hist->cutlleePassHistW->Fill(cutFlowCH::Trigger, weight);

	// Ensure that event has the minimum number of required leptons: 2e and 2mu
	if(!(elLooseEvent.size() >= 2 &&  muEvent.size() >= 2 ) && !(elLooseEvent.size() >= 4)) return false;
	cutlleePass[cutFlowCH::Lepton]++;
	Hist->cutlleePassHistW->Fill(cutFlowCH::Lepton, weight);
	
	// Getting the diLepton combination
	diEvent2El.clear();
	diEvent2El = GetDiLeptonComb(elLooseEvent);
	diEvent2Mu.clear();
	diEvent2Mu = GetDiLeptonComb(muEvent);
	if(!(diEvent2El.size() >= 1 && diEvent2Mu.size() >= 1) && !(diEvent2El.size() >= 2)) return false;

	// Get a single Quad Event based on SFOS and PGD mass
	vector<QuadLepton *> higgsContainer =  GetQuadLeptonComb(diEvent2Mu, diEvent2El /*, quadType::Mu2El2 */);
	vector<QuadLepton *> higgsContainer4e =  GetQuadLeptonComb(diEvent2El);
		
	vector<QuadLepton*> combinedHiggsContainer = higgsContainer;
	combinedHiggsContainer.insert(combinedHiggsContainer.end(), higgsContainer4e.begin(), higgsContainer4e.end());
	eventToCleanCR = combinedHiggsContainer;	
	
	// Find the Z1 dilepton which passes the loose cut and it is closest to Z1
	// and Doesn't have overlaping looserelaxed leptons
	if(isDebugCallllee)
	{
		cout<<"--------------------------------------"<<endl;
		cout<<"Num of candidate quads: "<<combinedHiggsContainer.size()<<endl;
		cout<<"--------------------------------------"<<endl;
	}
	vector<QuadLepton*> cleanZ14e;
	//cleanZ14e = cutClosestZ1(higgsContainer);
	vector<QuadLepton*> cleanZ12e2mu;
	//cleanZ12e2mu = cutClosestZ1(higgsContainer4e);
	
	vector<QuadLepton*> cleanZ1 = cleanZ14e;
	cleanZ1.insert(cleanZ1.end(), cleanZ12e2mu.begin(), cleanZ12e2mu.end());


	// to ensure that I don't pick 2mu2e over 4e
	cleanZ1.clear();
	cleanZ1 = cutClosestZ1(combinedHiggsContainer);
	
	//cout<<"diEvent2El: "<<diEvent2El.size()<<" diEvent2Mu: "<<diEvent2Mu.size()<<endl;
	//cout<<"Before clean: "<<combinedHiggsContainer.size()<<" after: "<<cleanZ1.size()<<endl<<endl;

	if(cleanZ1.size() == 0) return false;
	cutlleePass[cutFlowCH::SFOS]++;
	Hist->cutlleePassHistW->Fill(cutFlowCH::SFOS, weight);

	if(isDebugCallllee)
	{
		cout<<"--------------------------------------"<<endl;
		cout<<"Num of after clean quads: "<<cleanZ1.size()<<endl;
		cout<<"--------------------------------------"<<endl;
	}
	for(Int_t i = 0; i < (Int_t) cleanZ1.size(); i++)
	{
		QuadLepton *currHiggs = cleanZ1[i];

		// Cut Quad lepton on kinematics and mass
		// This is a hack to use the same functions on these events
		vector<QuadLepton*> vecHiggs;
		vecHiggs.push_back(currHiggs);
		
		//QuadLepton* higgs = GetQuadEvent(vecHiggs, cutlleePass, currHiggs->getQuadType(), Hist->cutlleePassHistW, weight);
		///if(higgs == 0) continue;
		
		QuadLepton* higgs = currHiggs;
		Bool_t passCut = CutQuadLepton(higgs, cutlleePass, currHiggs->getQuadType(), Hist->cutlleePassHistW, &weight,prodCHllee, true);
		vecHiggs.clear();
		if(!passCut) continue;
		
		// Saving the weight
		higgs->weight_corr = weight;
	
		outputTreelleeCR->fillTree(event, higgs,  higgs->getQuadType(), isMC);
		//fillEventVarInfo(higgsCandidate2L2L, type, prodCH2L2L);
	}
	if(isDebugCallllee) isDebugCall = false;

}

// To fill the information on the higgs
void HiggsAnalysis::fillEventVarInfo(QuadLepton * higgs, Int_t type, Int_t *prodCH)
{
	// FSR Correction
	MassCalc(higgs, type);
	if(isDebugCall | printWeight)
	{
		cout<<"----------------"<<endl;
		cout<<"Higgs Mass: "<< higgs->getMass()<<endl;
		cout<<"Higgs Mass err: "<< higgs->getMassErr()<<endl;
		cout<<"Higgs FSR Mass: "<<higgs->getMassFSR()<<endl;
		cout<<"Higgs FSR Mass err: "<<higgs->getMassErrFSR()<<endl;
		cout<<"Higgs Constrained Mass: "<<higgs->getMassZMassCons()<<endl;		
		cout<<"Higgs Constrained Mass err: "<<higgs->getMassErrZMassCons()<<endl;
		cout<<"m4l_unconstrained_ID 	"<<higgs->massID/1000<<endl;
		cout<<"m4lerr_unconstrained_ID 	"<<higgs->massErrID/1000<<endl;
		cout<<"m4l_constrained_ID 		"<<higgs->massZMassConsID/1000<<endl;
		cout<<"m4lerr_constrained_ID	"<<higgs->massErrZmassConsID/1000<<endl;
		cout<<"m4l_unconstrained_MS 	"<<higgs->massMS/1000<<endl;
		cout<<"m4lerr_unconstrained_MS 	"<<higgs->massErrMS/1000<<endl;
		cout<<"m4l_constrained_MS 		"<<higgs->massZMassConsMS/1000<<endl;
		cout<<"m4lerr_constrained_MS	"<<higgs->massErrZmassConsMS/1000<<endl;

		cout<<"----------------"<<endl;

	}
	// Get Production Channel
	fillProductionChannel(higgs, muEvent, elEvent, jetsEvent);
	
	// Angular Vars
	GetElicityAngles(higgs->getZ1()->getLepNeg()->get4MomentumNoP(),
					higgs->getZ1()->getLepPlus()->get4MomentumNoP(),
					higgs->getZ2()->getLepNeg()->get4MomentumNoP(),
					higgs->getZ2()->getLepPlus()->get4MomentumNoP(),
					&higgs->cthstr,
					&higgs->phi1,
					&higgs->cth1,
					&higgs->cth2,
					&higgs->phi);
	
	// Truth fill
	if(generatorName == MCGeneratorName::Pythia) fillTruthRecoMatchedInfo(higgs);

	// Filling weight
	fillHiggsWeight(higgs);
	fillCrossSection(higgs);
	
	// Flag for dataOverlap
	fillFlagStream(higgs);

	// LepId
	fillLepID(higgs);

	// Electron ID
	fillLooseCut(higgs);
	
	// Filling the true jets info
	fillTruthJetsInfo(higgs);

	// Fill jet infor for fudicial people
	fillFudicialJets(higgs);

	// BDT stuff
	fillBDTWeights(higgs);

	// NPV
	higgs->npv = getNVertex(2); 

	// Calibration type
	higgs->calib = curCalibration;

}
// Choose all the quadleptons that are closest to zPDG and are loose electrons
// plus removes event that are overlapping in etaphi electrons
vector<QuadLepton*> HiggsAnalysis::cutClosestZ1(vector<QuadLepton*> higgsContainer)
{
	vector<QuadLepton *> z1Overlap;

	for(Int_t i = 0; i < (Int_t) higgsContainer.size(); i++)
	{
		QuadLepton* currHiggs = higgsContainer[i];
		//cout<<"currHiggs: "<<currHiggs<<" i: "<<i<<" size: "<<higgsContainer.size()<<endlu;
		vector<ChargedLepton *> currLep = currHiggs->getLepton();

		//cout<<currLep.size()<<endl;
		Bool_t isOverlap = false;
		// For each lepton in container
		for(Int_t j = 0; j < (Int_t) currLep.size(); j++)
		{
			if(currLep[j]->getFlavor() == flavor::Muon) continue;

			// else loop over over the continer to find overlap
			for(Int_t k = 0; k <(Int_t) currLep.size(); k++)
			{
				if(currLep[k]->getFlavor() == flavor::Muon) continue;
				if(j == k) continue;
				
				// Track overlap
				D3PDReader::ElectronD3PDObjectElement *el_curr = currLep[j]->GetElectron();
				// Getting the Current Variables
				Double_t d0_curr = -1;
				Double_t z0_curr = -1;
				Double_t phi_curr = -1;
				Double_t qoverp_curr = -1;
				Double_t Eta_curr = -1;

				if(dataYear == 2011)
				{
					d0_curr = el_curr->trackd0();
					z0_curr = el_curr->trackz0();
					phi_curr = el_curr->trackphi();
					qoverp_curr = el_curr->trackqoverp();
					Eta_curr = el_curr->tracketa();
				}
				else if (dataYear == 2012)
				{
					d0_curr = el_curr->Unrefittedtrack_d0();
					z0_curr = el_curr->Unrefittedtrack_z0();
					phi_curr = el_curr->Unrefittedtrack_phi();
					qoverp_curr = el_curr->Unrefittedtrack_qoverp();
					Eta_curr = el_curr->Unrefittedtrack_eta();
				}
				D3PDReader::ElectronD3PDObjectElement *el_j = currLep[j]->GetElectron();
				if(el_curr == el_j) continue;

				/// For the electron to compare
				Double_t d0_j = -1;
				Double_t z0_j = -1;
				Double_t phi_j = -1;
				Double_t qoverp_j = -1;
				Double_t Eta_j = -1;

				if(dataYear == 2011)
				{
					d0_j = el_j->trackd0();
					z0_j = el_j->trackz0();
					phi_j = el_j->trackphi();
					qoverp_j = el_j->trackqoverp();
					Eta_j = el_j->tracketa();
				}
				else if (dataYear == 2012)
				{
					d0_j = el_j->Unrefittedtrack_d0();
					z0_j = el_j->Unrefittedtrack_z0();
					phi_j = el_j->Unrefittedtrack_phi();
					qoverp_j = el_j->Unrefittedtrack_qoverp();
					Eta_j = el_j->Unrefittedtrack_eta();
				}

				// Comparing they share the same ID and if the curr has lower ET, reject it
				if( d0_curr == d0_j &&
					z0_curr == z0_j &&
					phi_curr == phi_j &&
					qoverp_curr == qoverp_j
				  ) isOverlap = true;


				// Calo overlap
				Double_t phi_Clcurr = -1;
				Double_t eta_Clcurr = -1;

				Double_t cutPhi = -1;
				Double_t cutEta = -1;

				if(dataYear == 2011)
				{
					cutPhi = 0;
					cutEta = 0;
				}
				else if(dataYear == 2012)
				{
					cutPhi = 5*0.025;
					cutEta = 3*0.025;
				}

				phi_Clcurr = el_curr->cl_phi();
				eta_Clcurr = el_curr->cl_eta();

				Double_t phi_Clj = el_j->cl_phi();
				Double_t eta_Clj = el_j->cl_eta();

				Double_t deltaEta = fabs(eta_Clcurr - eta_Clj);
				Double_t deltaPhi = fabs(phi_Clcurr - phi_Clj);
				deltaPhi = (deltaPhi > TMath::Pi()) ? 2*TMath::Pi()-deltaPhi : deltaPhi;

				// Comparing they share the same ID and if the curr has lower ET, reject it
				if( deltaEta < cutEta &&
						deltaPhi < cutPhi 
				  ) isOverlap = true;;

			}
		}
		vector<QuadLepton*> vecHiggs;
		vecHiggs.push_back(higgsContainer[i]);
		Double_t weight = 1;	
		QuadLepton* higgs = GetQuadEvent(vecHiggs, cutlleePass, currHiggs->getQuadType(), Hist->cutlleePassHistW, weight);
		vecHiggs.clear();
		if(isDebugCallllee)
		{
			cout<<"--------------------------------------"<<endl;
			if(isOverlap) cout<<"Reject Quad at electron Overlap"<<endl;
			if(higgs == 0) cout<<"Reject Quad as is it didn't pass the GetQuadEvent cuts"<<endl;
			cout<<"--------------------------------------"<<endl;
		}

		if(higgs == 0) continue;
		if(!isOverlap) z1Overlap.push_back(higgsContainer[i]);
	}
	// LepID cut
	vector<QuadLepton*> z1IDCut;

	for(Int_t i = 0; i < (Int_t) z1Overlap.size(); i++)
	{
		if(isDebugCallllee) isDebugCall = false;

		// filling the ID 
		fillLooseCut(z1Overlap[i]);

		if(isDebugCallllee) isDebugCall = true;

		Bool_t idCut = z1Overlap[i]->looseElectron[0] & z1Overlap[i]->looseElectron[1];

		// if Z1 passes the ID, save the quad
		if(idCut) z1IDCut.push_back(z1Overlap[i]);
		else
		{
			if(isDebugCall) cout<<"Rejected at Z1 id Cut"<<endl;	
		}
	}

	// closes Z1 cut to ZPDG mass
	Double_t diffZmin = 9999999;
	DiLepton* minZ = NULL;
	for(Int_t i = 0; i < (Int_t) z1IDCut.size(); i++)
	{
		Double_t diffCurr = fabs(z1IDCut[i]->getZ1()->get4Momentum()->M() - pdgZMass);
		if(diffCurr < diffZmin)
		{
			diffZmin = diffCurr;
			minZ = z1IDCut[i]->getZ1();
		}
	}

	vector<QuadLepton *> z1MassCut;
	for(Int_t i = 0; i < (Int_t) z1IDCut.size(); i++)
	{
		if(minZ == z1IDCut[i]->getZ1() && z1IDCut[i]->getZ2()->getFlavor() != flavor::Muon) z1MassCut.push_back(z1IDCut[i]);
	}

	// Cut based on iso and d0 sig
	vector<QuadLepton*> finalZ1Quad;
	for(Int_t i = 0; i < (Int_t) z1MassCut.size(); i++)
	{
		//fillLooseCut(z1MassCut[i]);

		if(isDebugCallllee) isDebugCall = false;
		vector<Double_t> pTCone20 = getTrackPT(z1MassCut[i]->getLepton());
		Bool_t passCutTrackIso = CutTrackIso(z1MassCut[i]->getLepton(), pTCone20, z1MassCut[i]->getQuadType(), z1MassCut[i]);

		vector<Double_t> eTCone20 = getCaloET(z1MassCut[i]->getLepton());
		Bool_t passCutCaloIso = CutCaloIso(z1MassCut[i]->getLepton(), eTCone20, z1MassCut[i]->getQuadType(), z1MassCut[i]);

		Bool_t passCutD0Sig = D0SigCut(z1MassCut[i]->getLepton(), z1MassCut[i]->getQuadType(), z1MassCut[i]);
		if(isDebugCallllee) isDebugCall = true;

		Bool_t authorCut = false;
		authorCut 		= z1MassCut[i]->looseElectron[0] & z1MassCut[i]->looseElectron[1];
		passCutTrackIso = z1MassCut[i]->trackIso[0] & z1MassCut[i]->trackIso[1];
		passCutCaloIso 	= z1MassCut[i]->caloIso[0] & z1MassCut[i]->caloIso[1];
		passCutD0Sig	= z1MassCut[i]->d0Sig[0] & z1MassCut[i]->d0Sig[1];
		if(isDebugCallllee)
		{
			cout<<"---------------------------"<<endl;
			cout<<"llee Cutflow"<<endl;
			for(Int_t j = 0; j < z1MassCut[i]->looseElectron.size(); j++) cout<<"LepID "<<j <<" :"<< z1MassCut[i]->looseElectron[j]<<endl;
			for(Int_t j = 0; j < z1MassCut[i]->trackIso.size(); j++) cout<<"trackIso "<<j <<" :"<< z1MassCut[i]->trackIso[j]<<" "<<z1MassCut[i]->valTrackIso[j]<<endl;
			for(Int_t j = 0; j < z1MassCut[i]->caloIso.size(); j++) cout<<"caloIso "<<j <<" :"<< z1MassCut[i]->caloIso[j]<<" "<<z1MassCut[i]->valCaloIso[j]<<endl;
			for(Int_t j = 0; j < z1MassCut[i]->d0Sig.size(); j++) cout<<"d0Sig "<<j <<" :"<< z1MassCut[i]->d0Sig[j]<<" "<<z1MassCut[i]->valD0Sig[j]<<endl;
		}
		if(authorCut && passCutTrackIso && passCutCaloIso && passCutD0Sig) finalZ1Quad.push_back(z1MassCut[i]);
		else
		{
			if(isDebugCallllee)
			{
				cout<<"--------------------------------------"<<endl;
				cout<<"Reject Quad at Z1 trackiso and sig cut"<<endl;
				cout<<"--------------------------------------"<<endl;
			}
		}
	}
	z1MassCut.clear();
	z1IDCut.clear();

	return finalZ1Quad;

}

void HiggsAnalysis::fillLooseCut(QuadLepton* higgs)
{

	higgs->looseElectron.clear();
	vector<ChargedLepton *> curr_lep = higgs->getLepton();
	
	for(UInt_t i = 0; i < curr_lep.size(); i++)
	{
		ChargedLepton * lep_i = curr_lep[i];
		if(lep_i->getFlavor() == flavor::Electron)
		{
			D3PDReader::ElectronD3PDObjectElement *el_i = lep_i->GetElectron();

			higgs->looseElectron.push_back(electronObj->AuthorCut(el_i,lep_i->smearVal, getNVertex(2), false));
		}
		if(lep_i->getFlavor() == flavor::Muon)
		{
			higgs->looseElectron.push_back(true);
		}

	}

	if(higgs->looseElectron.size() != 4) cout<<"Error: size of looseElectron vector not 4"<<endl;
	curr_lep.clear();

}

////////////////////////////////////////////////////////////////////////////////////////
//			Initialize the variables and Selector Tools
////////////////////////////////////////////////////////////////////////////////////////
// Initializes the variable for the each event
void HiggsAnalysis::InitializeEventVar()
{	
	// Checking if MC or Data
	if(event->eventinfo.isSimulation()) isMC = 1;
	else isMC = 0;

	// Finding the year and CM Energy of the data
	if((event->eventinfo.RunNumber()>=177531) && (event->eventinfo.RunNumber()<=191933))
	{
 		dataYear = 2011;
   		CM_E = 7.0;
		//curMCCollection = MCCollection::MC11c;
  	}
	else if (event->eventinfo.RunNumber()>191933 )
	{
   		dataYear = 2012;
   		CM_E = 8.0;
		//if(event->eventinfo.RunNumber() == 195847) curMCCollection = MCCollection::MC12a;
		//else if(event->eventinfo.RunNumber() == 195848) curMCCollection = MCCollection::MC12b;
		//else if(!isMC) curMCCollection = MCCollection::MC12a;
		//else cout<<"InitializeEventVar: currMCCollection not recognized"<<endl;
  	}
	else
	{
		cout<<"Event number not recognized: "<<event->eventinfo.RunNumber()<<endl;
	}

	// 2012 has GSF by deflaut. This is to take care of this
	if(dataYear == 2011) {el_cur = &(event->el_GSF);}
	else if(dataYear == 2012) { el_cur = &(event->el);}

	// Trigger Matching
	if(dataYear == 2011) D3DPTriggerDev = true;
	else if(dataYear == 2012) D3DPTriggerDev = false;

	// Geting the currMCcolliection
	FillProductionTag();
}


// Initialize all the pileup tool for analysis
void HiggsAnalysis::InitPileupTool()
{
	pileupTool = new Root::TPileupReweighting( "PileupReweightingTool" );

	event->GetEntry(0);
	InitializeEventVar();

	if(dataYear == 2011)
	{
		string name_PU="InputFile/PileupConfigFile/MC11c.prw.root";
		pileupTool->SetUnrepresentedDataAction(2);
		pileupTool->AddConfigFile(name_PU);
		pileupTool->AddLumiCalcFile("InputFile/PileupConfigFile/ilumicalc_2011_AllYear_All_Good.root"); 
		//pileupTool->AddLumiCalcFile("InputFile/PileupConfigFile/ilumicalc_2011_AllYear_All_Good.root_old"); 		
		pileupTool->SetDefaultChannel(109292);
		pileupTool->UsePeriodConfig("MC11c");
		pileupTool->Initialize();		
	}
	else if(dataYear == 2012 && (curMCCollection == MCCollection::MC12a || curDataCollection != -1))
	{
		string name_PU = "";
		if(doMoriondConfig) name_PU = "InputFile/PileupConfigFile/MC12a_Moriond.prw.root";
		else name_PU = "InputFile/PileupConfigFile/MC12a.prw.root";
		cout<<"Pileup Config file: "<<name_PU<<endl;
		pileupTool->SetUnrepresentedDataAction(2);
		pileupTool->AddConfigFile(name_PU);

		string name_Lumi = "";		
		if(doMoriondConfig) name_Lumi = "InputFile/PileupConfigFile/ilumicalc_2012_AllYear_All_Good_Moriond.root";
		else name_Lumi = "InputFile/PileupConfigFile/ilumicalc_2012_AllYear_All_Good.root";
		cout<<"Pileup lumi file: "<<name_Lumi<<endl;

		pileupTool->AddLumiCalcFile(name_Lumi); 
		pileupTool->SetDefaultChannel(160156);
		//pileupTool->AddBinning("pileup",50,0,50);
		pileupTool->Initialize();
	}
	else if(dataYear == 2012 && curMCCollection == MCCollection::MC12b)
	{
		string name_PU = "InputFile/PileupConfigFile/MC12b.prw.root";
		cout<<"Pileup Config file: "<<name_PU<<endl;
		
		pileupTool->SetUnrepresentedDataAction(2);
		pileupTool->AddConfigFile(name_PU);

		string name_Lumi = "InputFile/PileupConfigFile/ilumicalc_2012_AllYear_All_Good.root";
		cout<<"Pileup lumi file: "<<name_Lumi<<endl;

		pileupTool->AddLumiCalcFile(name_Lumi); 
		pileupTool->SetDefaultChannel(181341);
		//pileupTool->AddBinning("pileup",50,0,50);
		pileupTool->Initialize();
	}
	else if(dataYear == 2012 && curMCCollection == MCCollection::MC12c)
	{
		string name_PU = "InputFile/PileupConfigFile/MC12b.prw.root";
		cout<<"WARNING: NEED TO FIX THIS: USING MC12B SETUP FOR MC12c"<<endl;
		cout<<"Pileup Config file: "<<name_PU<<endl;
		
		pileupTool->SetUnrepresentedDataAction(2);
		pileupTool->AddConfigFile(name_PU);

		string name_Lumi = "InputFile/PileupConfigFile/ilumicalc_2012_AllYear_All_Good.root";
		cout<<"Pileup lumi file: "<<name_Lumi<<endl;

		pileupTool->AddLumiCalcFile(name_Lumi); 
		pileupTool->SetDefaultChannel(181341);
		//pileupTool->AddBinning("pileup",50,0,50);
		pileupTool->Initialize();
	}


	else cout<<"Error: InitPileUpTool: dataYear not recognized"<<endl;
}
void HiggsAnalysis::InitZVertexTool()
{
	if(dataYear == 2011)
	{
		if(event->eventinfo.mc_channel_number() == 105200)
		{
        	zVertexTool= new VertexPositionReweightingTool(std::string("s1272"),
					"../../egammaAnalysisUtils/share/zvtx_weights_2011_2012.root");		
		}
		else
		{
        	zVertexTool= new VertexPositionReweightingTool(std::string("s1310"),
					"../../egammaAnalysisUtils/share/zvtx_weights_2011_2012.root");
		}
	}
	else if(dataYear == 2012)
	{
   		zVertexTool= new VertexPositionReweightingTool(VertexPositionReweightingTool::MC12a,
				"../../egammaAnalysisUtils/share/zvtx_weights_2011_2012.root");
	}
}
// Good Run list
void HiggsAnalysis::InitGoodRunList()
{
	TString goodrunslistname;

	if(dataYear == 2011) goodrunslistname = "InputFile/GRL/data11_7TeV.periodAllYear_DetStatus-v36-pro10-02_CoolRunQuery-00-04-08_All_Good.xml";
	else if(dataYear == 2012)
	{
		if(doMoriondConfig) goodrunslistname = "InputFile/GRL/data12_8TeV.periodAllYear_DetStatus-v58-pro14-01_DQDefects-00-00-33_PHYS_StandardGRL_All_Good.xml";
		else goodrunslistname = "InputFile/GRL/data12_8TeV.periodAllYear_DetStatus-v61-pro14-02_DQDefects-00-01-00_PHYS_StandardGRL_All_Good.xml";
		//goodrunslistname = "../../GRL/data12_8TeV.periodAllYear_DetStatus-v61-pro14-02_DQDefects-00-01-00_PHYS_StandardGRL_All_Good.xml";
	}
	else {cout<<"Error: HiggsAnalysis::InitGoodRunList: dataYear not recognized"<<endl;}
	
	Root::TGoodRunsListReader* GoodRunList = new Root::TGoodRunsListReader();
	cout<<"GRL list: "<< goodrunslistname<<endl;
    GoodRunList->AddXMLFile(goodrunslistname);
	GoodRunList->Interpret();
	grl = GoodRunList->GetMergedGoodRunsList();
	
	
}
// TileTrip
void HiggsAnalysis::InitTileTrip()
{
	TileTrip=new Root::TTileTripReader("TripReader");
	TileTrip->setTripFile("../../TileTripReader/data/CompleteTripList_2011-2012.root>" );
}

// ggF Reweighing
void HiggsAnalysis::InitggFReweight()
{
	cout<<"Going for ggF reweight"<<endl;
	// Don't init if it data as it requires sample mass and the fact that sample is ggF
	if(!isMC || !doggFWeight)
	{

		ggFReweight = 0;
		cout<<"ggF would have exited here: isMC: "<<isMC<<" doggFWeight: "<<doggFWeight<<endl;
		//return;
	}
	Double_t samplemass = getMCHiggsMass();
	// Don't need for not ggF samples
	if(sampleProdType != sampleType::ggF && sampleProdType != sampleType::ggF_ZpZp)
	{
		ggFReweight = 0;
		cout<<"ggF exited here 2: not found to be a ggH sample"<<endl;	
		return;
	}
	Double_t massHiggsggFR = samplemass;
	Double_t DeltaMggFR = 9999;

	Int_t powheg_mass[] =  {100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155,
                                160, 165, 170, 175, 180, 185, 190, 195, 200, 210, 220, 240,
                                260, 280, 300, 320, 340, 360, 380, 400, 420, 440, 460, 480,
                                500, 520, 540, 560, 580, 600, 650, 700, 750, 800, 850, 900 };
 
	for (Int_t i = 0; i<(Int_t) (sizeof((powheg_mass))/sizeof((powheg_mass[0])));i++) {
		Double_t tmpDeltaM = double(powheg_mass[i])-samplemass;
		if ( fabs(tmpDeltaM)<fabs(DeltaMggFR) ){
			//std::cout <<" M tool -M sample = "<< tmpDeltaM << " Mass tool = " <<  double(powheg_mass[i]) <<endl;
			DeltaMggFR    = tmpDeltaM;
			massHiggsggFR = powheg_mass[i];
		}
	}
	cout<<"Setting up ggFReweight, sample mass: "<<massHiggsggFR<<endl;
	ggFReweight = new ggFReweighting("PowHeg", massHiggsggFR, "Mean", "../../ggFReweighting/share/", "mc11");
	
}
// BDT tool for categorizations
void HiggsAnalysis::InitCategoryBDTTool()
{

	dijet_invmass = -999;
	dijet_deltaeta = -999;
	leading_jet_pt = -999;
	leading_jet_eta = -999;
	subleading_jet_pt = -999;

	BDT_discriminant_VBF = -999;
	BDT_discriminant_HadVH = -999;

	cout<<"Initialzing the CategoriesDiscriminantTool"<<endl<<endl;
	CategoriesDiscriminantTool = new CategoriesMVA("../../CategoriesMVA/weights");

	CategoriesDiscriminantTool->ConnectVariables  (dijet_invmass, dijet_deltaeta, leading_jet_pt, leading_jet_eta, subleading_jet_pt);
	
}

void HiggsAnalysis::InitBDTTool()
{
	BDTtool = new H4lBDTWeights();
	//cout.setf(ios::fixed, ios::floatfield);
	//cout.setf(ios::showpoint);
}

// Pt reweighing for JHU samples
void HiggsAnalysis::InitJHUReweight()
{
	if(dataYear == 2011)
	{
      ptJHUReweighting= new JHUPtReweighting("mc11","../../JHUReweighting/share/");
	}
	else if(dataYear == 2012)
	{
      ptJHUReweighting= new JHUPtReweighting("mc12","../../JHUReweighting/share/");		
	}
	else {cout<<"InitJHUReweigh: DataYear not recognized"<<endl;}
}


// Helper for initializing
void HiggsAnalysis::getPeriodEvent()
{
	Int_t nRunTemp = event->eventinfo.RunNumber();
	//2011
	//// 2011 period B-D
  	if(nRunTemp>=177986 && nRunTemp<=180481) dataPeriod = period2011_BD; 
  	//// 2011 period E-F-G-H
  	if(nRunTemp>=180614 && nRunTemp<=184169) dataPeriod = period2011_EH;
	//// 2011 period I-J-K
  	if(nRunTemp>=185353 && nRunTemp<=187815) dataPeriod = period2011_IK;
  	//// 2011 period L
  	if(nRunTemp>=188902 && nRunTemp<=191933) dataPeriod = period2011_LM;
 	//////// 2012
  	//// 2012 period All 
  	if(nRunTemp>=195847 && nRunTemp<=999999) dataPeriod = period2012_All;

	// period I [185353-186493]
	//if(nRunTemp>=185353 && nRunTemp<=186493) { nI++;}
	// period J [186516-186755]
     	//if(nRunTemp>=186516 && nRunTemp<=186755) { nJ++;}
	// period K [186873-187815
     	//if(nRunTemp>=186873 && nRunTemp<=187815) { nK++;}
	//std::cout<<dataPeriod<<endl;
}
////////////////////////////////////////////////////////////////////////////////////////
//				CutFlow Functions
//		all Function return 1 or true if event passed the cut
////////////////////////////////////////////////////////////////////////////////////////
Bool_t HiggsAnalysis::DataPreselectionCut()
{	
	if(event->eventinfo.larError() == 2) return false;
	if(event->eventinfo.tileError() == 2) return false;
	if((event->eventinfo.coreFlags() & 0x40000) != 0) return false;
	if(! (TileTrip->checkEvent(event->eventinfo.RunNumber(), event->eventinfo.lbn(), event->eventinfo.EventNumber()))) return false;
	// Good Run List
	int ReinitialzedTool=grl.HasRunLumiBlock(event->eventinfo.RunNumber(),-1); 		
	if(grl.HasRunLumiBlock(event->eventinfo.RunNumber(),event->eventinfo.lbn()))
	{}
	else { return false; }


	return true;
}

Bool_t HiggsAnalysis::AllPreselectionCut()
{
	Bool_t passCut = true;
	passCut = VertexCut();
	if(!passCut) return passCut;

	return passCut;
}
////////////////////////////////////////////////////////////////////////////////////////
//				Helper CutFlow Functions
//		all Function return 1 or true if event passed the cut
////////////////////////////////////////////////////////////////////////////////////////
// Cuts on the fact that atleast one vertex must have more than nMinVertes tracks
Bool_t HiggsAnalysis::VertexCut()
{
	int nVertex = event->vxp.n();
	for (int i = 0; i < nVertex; i++)
	{
		if(event->vxp[i].trk_n() >= 3) return true;
	}
	return false;
}
Int_t HiggsAnalysis::getNVertex(Int_t nCutVertex)
{
	Int_t nVertex = 0; 
	for (int i = 0; i < event->vxp.n(); i++)
	{
		if(event->vxp[i].trk_n() >= nCutVertex) nVertex++;
	}
	return nVertex;
}

// Single Electron Trigger
Bool_t HiggsAnalysis::SingleElectronTrigger(Int_t runNumber)
{
	Bool_t pass = true;
	if(dataPeriod == period2011_BD || dataPeriod == period2011_EH) pass = event->triggerbits.EF_e20_medium();
	
	// For the random period...
	else if(dataPeriod == period2011_IK)
	{
		if(runNumber < 186873) pass = event->triggerbits.EF_e20_medium();
		else pass = event->triggerbits.EF_e22_medium();
	}
	else if(dataPeriod == period2011_LM) pass = event->triggerbits.EF_e22vh_medium1();
	else if(dataPeriod == period2012_All) pass = event->triggerbits.EF_e24vhi_medium1() |
						event->triggerbits.EF_e60_medium1();
	return pass;
}

// Di Electron Trigger
Bool_t HiggsAnalysis::DiElectronTrigger(Int_t runNumber)
{
	Bool_t pass = true;
	if(dataPeriod == period2011_BD || dataPeriod == period2011_EH) pass = event->triggerbits.EF_2e12_medium();
	else if(dataPeriod == period2011_LM) pass = event->triggerbits.EF_2e12Tvh_medium();
	else if(dataPeriod == period2012_All) 
	{
		if(!isMC) pass = event->triggerbits.EF_2e12Tvh_loose1() | event->triggerbits.EF_2e12Tvh_loose1_L2StarB();
		else pass = event->triggerbits.EF_2e12Tvh_loose1();
	}

	// random period...	
	else if(dataPeriod == period2011_IK)
	{
		if(runNumber < 186873) pass = event->triggerbits.EF_2e12_medium();
		else pass = event->triggerbits.EF_2e12T_medium();
	}


	return pass;
}
// Single Muon Trigger
Bool_t HiggsAnalysis::SingleMuonTrigger(Int_t runNumber)
{
	Bool_t pass = true;
	if(dataPeriod == period2011_BD || dataPeriod == period2011_EH 
	) pass = event->triggerbits.EF_mu18_MG();
	else if(dataPeriod == period2011_LM) pass = event->triggerbits.EF_mu18_MG_medium();
	else if(dataPeriod == period2012_All) pass = event->triggerbits.EF_mu24i_tight() |
					 event->triggerbits.EF_mu36_tight();
	// random period...
	else if(dataPeriod == period2011_IK)
	{
		if(runNumber < 186516 ) pass = event->triggerbits.EF_mu18_MG();
		else  pass = event->triggerbits.EF_mu18_MG_medium();
	}				 		
	return pass;
}
// Di Muon Trigger
Bool_t HiggsAnalysis::DiMuonTrigger()
{
	Bool_t pass = true;
	if(dataPeriod == period2011_BD || dataPeriod == period2011_EH 
	|| dataPeriod == period2011_IK 
	|| dataPeriod == period2011_LM) pass = event->triggerbits.EF_2mu10_loose();
	if(dataPeriod == period2012_All) pass = event->triggerbits.EF_2mu13() |
					 	event->triggerbits.EF_mu18_tight_mu8_EFFS();
	return pass;

}
//E-Mu trigger
Bool_t HiggsAnalysis::ElectronMuonTrigger()
{
	Bool_t pass = true;
	if(dataPeriod == period2011_BD || dataPeriod == period2011_EH 
	|| dataPeriod == period2011_IK 
	|| dataPeriod == period2011_LM) pass = event->triggerbits.EF_e10_medium_mu6();
	if(dataPeriod == period2012_All) pass = event->triggerbits.EF_e12Tvh_medium1_mu8() |
					 	event->triggerbits.EF_e24vhi_loose1_mu8();
	return pass;

}
// To remover e-mu, calo-e and jets-e overlap.
void HiggsAnalysis::RemoveOverlap()
{
	// El mu overlap
	for(vector<ChargedLepton *>::iterator itr = elOverlap.begin();
			itr != elOverlap.end(); ++itr)
	{
		ChargedLepton *el_lep_i = *itr;
		D3PDReader::ElectronD3PDObjectElement *el_i = el_lep_i->GetElectron();

		Bool_t reject = false;

		reject = ElectronMuonOverlap(el_i);

		if(!reject)
		{
			elEvent.push_back(el_lep_i);
			cutElPass[cutElFlow::OverLap] ++;
		}
	}
	// El mu overlap for loose electron
	for(vector<ChargedLepton *>::iterator itr = elLooseOverlap.begin();
			itr != elLooseOverlap.end(); ++itr)
	{
		ChargedLepton *el_lep_i = *itr;
		D3PDReader::ElectronD3PDObjectElement *el_i = el_lep_i->GetElectron();

		Bool_t reject = false;

		reject = ElectronMuonOverlap(el_i);
		
		if(! electronObj->AuthorCut(el_i,el_lep_i->smearVal, getNVertex(2), false)) reject = false;

		if(!reject)
		{
			elLooseEvent.push_back(el_lep_i);
			cutElLoosePass[cutElFlow::OverLap] ++;
		}
	}
	// Calo muon and el overlap
	for(vector<ChargedLepton *>::iterator itr = muOverlap.begin();
			itr != muOverlap.end(); ++itr)
	{
		ChargedLepton *mu_lep_i = *itr;
		D3PDReader::MuonD3PDObjectElement *mu_i = mu_lep_i->GetMuon();

		Bool_t reject = false;
		
		if(mu_lep_i->type == leptonType::MuonCalo)
		{reject = CaloMuonElectronOverlap(mu_i);}

		if(!reject)
		{
			muEvent.push_back(mu_lep_i);
			cutMuPass[cutMuFlow::OverLap] ++;
		}
	}
	// jets and electrom
	for(vector<ChargedLepton *>::iterator itr = jetsOverlap.begin();
			itr != jetsOverlap.end(); ++itr)
	{
		ChargedLepton *jets_lep_i = *itr;
		D3PDReader::JetD3PDObjectElement *jets_i = jets_lep_i->GetJets();

		Bool_t reject = false;
		
		reject = JetsElectronOverlap(jets_i);

		if(!reject)
		{
			jetsEvent.push_back(jets_lep_i);
			cutJetsPass[cutJetsFlow::OverLap] ++;
		}
	}

	// Fudicial jets and electrom
	for(vector<ChargedLepton *>::iterator itr = jetsOverlap_Fid.begin();
			itr != jetsOverlap_Fid.end(); ++itr)
	{
		ChargedLepton *jets_lep_i = *itr;
		D3PDReader::JetD3PDObjectElement *jets_i = jets_lep_i->GetJets();

		Bool_t reject = false;
		
		reject = JetsElectronOverlap(jets_i);

		if(!reject)
		{
			jetsEvent_Fid.push_back(jets_lep_i);
		}
	}

}
// Event Comparsion for electron overlap with muon
Bool_t HiggsAnalysis::ElectronMuonOverlap(D3PDReader::ElectronD3PDObjectElement *el_curr)
{
	Double_t el_eta = -1;
	Double_t el_phi = -1;

	if(dataYear == 2011)
	{
		el_eta = el_curr->tracketa();
		el_phi = el_curr->trackphi();
	}
	else if(dataYear == 2012)
	{
		el_eta = el_curr->Unrefittedtrack_eta();
		el_phi = el_curr->Unrefittedtrack_phi();
	}
	// Looping over to find it there is an overlap muon
	for(vector<ChargedLepton *>::iterator itr = muonObj->muBfOverlap.begin();
			itr != muonObj->muBfOverlap.end(); ++itr)
	{
		ChargedLepton *mu_lep_i = *itr;
		D3PDReader::MuonD3PDObjectElement *mu_i = mu_lep_i->GetMuon();
		if(mu_lep_i->type != leptonType::MuonCalo)
		{
			Double_t mu_theta = mu_i->id_theta();
			Double_t mu_eta = -log(tan(mu_theta*0.5));
			Double_t mu_phi = mu_i->id_phi();
			if(DeltaR(mu_eta, mu_phi, el_eta, el_phi) < 0.02) {return true;}
		}
	}

	return false;
}
// Overlap between Calo and electron
Bool_t HiggsAnalysis::CaloMuonElectronOverlap(D3PDReader::MuonD3PDObjectElement *mu_curr)
{
	Double_t mu_theta = mu_curr->id_theta();
	Double_t mu_eta = -log(tan(mu_theta*0.5));
	Double_t mu_phi = mu_curr->id_phi();

	// Looping over to find it there is an overlap electron
	for(vector<ChargedLepton *>::iterator itr = electronObj->elBfOverlap.begin();
			itr != electronObj->elBfOverlap.end(); ++itr)
	{
		ChargedLepton *el_lep_i = *itr;
		D3PDReader::ElectronD3PDObjectElement *el_i = el_lep_i->GetElectron();
		Double_t el_eta = -1;
		Double_t el_phi = -1;

		if(dataYear == 2011)
		{
			el_eta = el_i->tracketa();
			el_phi = el_i->trackphi();
		}
		else if(dataYear == 2012)
		{
			el_eta = el_i->Unrefittedtrack_eta();
			el_phi = el_i->Unrefittedtrack_phi();
		}
		if(DeltaR(mu_eta, mu_phi, el_eta, el_phi) < 0.02) {return true;}
	}

	return false;
}
// Overlap between jets and electron
Bool_t HiggsAnalysis::JetsElectronOverlap(D3PDReader::JetD3PDObjectElement *jets_curr)
{
	Double_t jets_eta = jets_curr->emscale_eta();
	Double_t jets_phi = jets_curr->emscale_phi();

	// Looping over to find it there is an overlap electron 
	for(vector<ChargedLepton *>::iterator itr = elEvent.begin();
			itr != elEvent.end(); ++itr)
	{
		ChargedLepton *el_lep_i = *itr;
		D3PDReader::ElectronD3PDObjectElement *el_i = el_lep_i->GetElectron();
		Double_t el_eta = -1;
		Double_t el_phi = -1;

		el_eta = el_i->tracketa();
		el_phi = el_i->trackphi();
		if(DeltaR(jets_eta, jets_phi, el_eta, el_phi) < 0.2) {return true;}
	}

	return false;
}
// This function create Dilepton pairs for the given set of leptons
vector<DiLepton *> HiggsAnalysis::GetDiLeptonComb(vector<ChargedLepton *>singleEvent)
{
	vector<DiLepton *> tempContainer;
	vector<DiLepton *> finalContainer;

	for(vector<ChargedLepton *>::iterator itr = singleEvent.begin();
			itr != singleEvent.end(); ++itr)
	{
		ChargedLepton *lep_i = *itr;
		// Looping over it again to find pairs
		for(vector<ChargedLepton *>::iterator itr_j = itr+1; 
			itr_j != singleEvent.end(); ++itr_j)
		{
			ChargedLepton *lep_j = *itr_j;
			DiLepton *temp;
			if(lep_j == lep_i) continue; // Don't want to compare the same thing

			// Have to be opposite charge
			if(! (lep_i->getCharge() * lep_j->getCharge() < 0)) continue;

			// The Dilepton class automatically orders the  
			temp = new DiLepton(lep_i, lep_j); 
			tempContainer.push_back(temp);
		}
	}
	if(isDebugCall) 
	{
		cout<<"-------------------"<<endl;
		cout<<"Number of Dilepton Candidate: "<<tempContainer.size()<<endl;
		cout<<"-------------------"<<endl;
	}
	
	return tempContainer;
}

// Returns the Quadlepton that fits the SFOS Criteria... for 4Mu and 4El
vector<QuadLepton *> HiggsAnalysis::GetQuadLeptonComb(vector<DiLepton *> diLepton)
{
	vector<QuadLepton *> higgsContainer;
	// Finding the primary lepton pair
	for(vector<DiLepton *>::iterator itr = diLepton.begin();
			itr != diLepton.end(); ++itr)
	{
		DiLepton *dilep_i = *itr;

		Int_t CaloCount = 0;
		Int_t StandAloneCount = 0;
		// Counting the calo and standalone muons
		if(dilep_i->getLepPlus()->getType() == leptonType::MuonCalo) {CaloCount++;}
		else if(dilep_i->getLepPlus()->getType() == leptonType::MuonStandAlone) {StandAloneCount++;}
		
		if(dilep_i->getLepNeg()->getType() == leptonType::MuonCalo) {CaloCount++;}
		else if(dilep_i->getLepNeg()->getType() == leptonType::MuonStandAlone) {StandAloneCount++;}

		for(vector<DiLepton *>::iterator itr_j = itr + 1;
			itr_j != diLepton.end(); ++itr_j)
		{
			DiLepton *dilep_j = *itr_j;
			if(dilep_j == dilep_i) continue;
			if(dilep_i->IsOverlap(dilep_j)) continue;

			Int_t CaloCountInt = 0;
			Int_t StandAloneCountInt = 0;

			if(dilep_j->getLepPlus()->getType() == leptonType::MuonCalo) {CaloCountInt++;}
			else if(dilep_j->getLepPlus()->getType() == leptonType::MuonStandAlone) {StandAloneCountInt++;}
			
			if(dilep_j->getLepNeg()->getType() == leptonType::MuonCalo) {CaloCountInt++;}
			else if(dilep_j->getLepNeg()->getType() == leptonType::MuonStandAlone) {StandAloneCountInt++;}

			// closest to z mass
			Double_t tempDiff_i = fabs(dilep_i->get4Momentum()->M() - pdgZMass);
			Double_t tempDiff_j = fabs(dilep_j->get4Momentum()->M() - pdgZMass);			
			
			if((CaloCount+StandAloneCount + CaloCountInt +StandAloneCountInt ) > 1) continue;	
				
			QuadLepton * temp = 0;	
			if(tempDiff_i < tempDiff_j)
				temp = new QuadLepton(dilep_i, dilep_j);
			else
				temp = new QuadLepton(dilep_j, dilep_i);

			higgsContainer.push_back(temp);
		}
	}

	if(isDebugCall) 
	{
		cout<<"-------------------"<<endl;
		cout<<"Number of 4l Quadlepton Candidate: "<<higgsContainer.size()<<endl;
		cout<<"-------------------"<<endl;
	}

	// To Find the secondary one
	return higgsContainer;

}

// Returns the Quadlepton that fits the SFOS Criteria... for 2Mu2El and 2El2Mu
vector<QuadLepton *> HiggsAnalysis::GetQuadLeptonComb(vector<DiLepton *> diMuLepton, vector<DiLepton *> diElLepton, Int_t forceQuadType)
{
	vector<QuadLepton *> higgsContainer;
	//cout<<"muContainer: "<<diMuLepton.size()<<" elContainer: "<< diElLepton.size()<<endl;
	// Finding the primary lepton pair
	for(vector<DiLepton *>::iterator itr = diMuLepton.begin();
			itr != diMuLepton.end(); ++itr)
	{
		DiLepton *dilep_i = *itr;

		Int_t CaloCount = 0;
		Int_t StandAloneCount = 0;
		// Counting the calo and standalone muons
		if(dilep_i->getLepPlus()->getType() == leptonType::MuonCalo) {CaloCount++;}
		else if(dilep_i->getLepPlus()->getType() == leptonType::MuonStandAlone) {StandAloneCount++;}
		
		if(dilep_i->getLepNeg()->getType() == leptonType::MuonCalo) {CaloCount++;}
		else if(dilep_i->getLepNeg()->getType() == leptonType::MuonStandAlone) {StandAloneCount++;}

		for(vector<DiLepton *>::iterator itr_j = diElLepton.begin();
			itr_j != diElLepton.end(); ++itr_j)
		{
			DiLepton *dilep_j = *itr_j;

			Int_t CaloCountInt = 0;
			Int_t StandAloneCountInt = 0; 

			if(dilep_j == dilep_i) continue;
			if(dilep_i->IsOverlap(dilep_j)) continue;
			
			if(dilep_j->getLepPlus()->getType() == leptonType::MuonCalo) {CaloCountInt++;}
			else if(dilep_j->getLepPlus()->getType() == leptonType::MuonStandAlone) {StandAloneCountInt++;}
			
			if(dilep_j->getLepNeg()->getType() == leptonType::MuonCalo) {CaloCountInt++;}
			else if(dilep_j->getLepNeg()->getType() == leptonType::MuonStandAlone) {StandAloneCountInt++;}
						
			// closest to z mass
			Double_t tempDiff_i = fabs(dilep_i->get4Momentum()->M() - pdgZMass);
			Double_t tempDiff_j = fabs(dilep_j->get4Momentum()->M() - pdgZMass);

			if((CaloCount+StandAloneCount + CaloCountInt +StandAloneCountInt ) > 1 && isDebugCall) cout<<"Reject as it exceed the calo and StandAlone count"<<endl ;	
			if((CaloCount+StandAloneCount + CaloCountInt +StandAloneCountInt ) > 1) continue;	
			
			QuadLepton * temp = 0;	
			if(tempDiff_i < tempDiff_j)
				temp = new QuadLepton(dilep_i, dilep_j);
			else
				temp = new QuadLepton(dilep_j, dilep_i);

			//cout<<"Quad m4l: "<<temp->get4Momentum()->M()<<" Z1: "<<temp->getZ1()->get4Momentum()->M();
			//cout<<" Z2: "<<temp->getZ2()->get4Momentum()->M();
			//if(temp->getQuadType() == analysisType::Mu2El2) cout<<" 2mu2el"<<endl;
			//else if(temp->getQuadType() == analysisType::El2Mu2) cout<<" 2el2mu"<<endl;
			//
			//vector<TLorentzVector> lepVec = temp->getLeptonLorentz();
			//for(Int_t i = 0; i < lepVec.size(); i++)
			//{
			//	cout<<"Lep: "<<i<<" pT: "<<lepVec[i].Pt()<<endl;
			//}
			// If we are forcing the quad type, not force it. 
			if(forceQuadType != -1)
			{
				if(temp->getQuadType() != forceQuadType) continue;
			}

			higgsContainer.push_back(temp);
		}

	}
	if(isDebugCall) 
	{
		cout<<"-------------------"<<endl;
		cout<<"Number of 2l2l Quadlepton Candidate: "<<higgsContainer.size()<<endl;
		cout<<"-------------------"<<endl;
	}

	return higgsContainer;

}
// Returns a single Quad Event based on SFOS And Trigger matching and PDG mass
QuadLepton*  HiggsAnalysis::GetQuadEvent(vector<QuadLepton *> higgs, Int_t * cutEventPass, Int_t type, TH1D *cutPassHistW, Double_t weight)
{

	QuadLepton *temp = 0;
	Double_t diffZ1Mass = 99999999;
	Double_t diffZ2Mass = 99999999;

  // Cut Values
	Double_t pTCut1 = 20*1000;
	Double_t pTCut2 = 15*1000;
	Double_t pTCut3 = 10*1000;

	Bool_t passKin = false;
	Bool_t passTrig = false;
	if(isDebugCall) 
	{
		cout<<"-------------------"<<endl;
		cout<<"Number of Higgs Candidate: "<<higgs.size()<<endl;
		cout<<"-------------------"<<endl;
	}
  	for(vector<QuadLepton *>::iterator itr = higgs.begin();
			itr != higgs.end(); ++itr)
	{
		int pTCount1 = 0;
		int pTCount2 = 0;
		int pTCount3 = 0;
		if(isDebugCall) cout<<"-------------------"<<endl<<"Quad Lepton Combinations Information: ";
		Double_t pTToStore [4] = {0};	

		QuadLepton* higgs_i = *itr;

		// Change the type depending on the current Higgs Candidate;
		if(higgs_i->getQuadType() == quadType::Mu2El2) {type = analysisType::Mu2El2;}
		else if(higgs_i->getQuadType() == quadType::El2Mu2) {type = analysisType::El2Mu2;}

		Int_t i = 0; // simple Counter
		// Loop over all the leptons inside the higgs class to find the 
		vector<ChargedLepton *> lep_i = higgs_i->getLepton();
		for(vector<ChargedLepton *>::iterator itr_lep = lep_i.begin();
			itr_lep != lep_i.end(); ++itr_lep)
		{
			ChargedLepton* curr_i = *itr_lep;
			Double_t pTCurr = curr_i->get4Momentum()->Pt();
			pTToStore[i] = pTCurr;
			i++;
			
			// For Z->4l analysis
			Double_t pTCut3Temp = pTCut3;
			if(anaType == doAnalysis::Zllll && curr_i->getFlavor() == flavor::Electron)  pTCut3Temp = 10 * 1000;
			else if(anaType == doAnalysis::Zllll && curr_i->getFlavor() == flavor::Muon) pTCut3Temp = 8 * 1000;

			if(pTCurr > pTCut1) pTCount1++;
			if(pTCurr > pTCut2) pTCount2++;
			if(pTCurr > pTCut3Temp) pTCount3++;
		}
		lep_i.clear();
		// To store the Pt
		std::sort(pTToStore, pTToStore + 4);
		Hist->quadInitPtHist[type][0]->Fill(pTToStore[0], Hist->weight);
		Hist->quadInitPtHist[type][1]->Fill(pTToStore[1], Hist->weight);
		Hist->quadInitPtHist[type][2]->Fill(pTToStore[2], Hist->weight);
		Hist->quadInitPtHist[type][3]->Fill(pTToStore[3], Hist->weight);
		// Cut
		if(pTCount1 >= 1 && pTCount2 >= 2 && pTCount3 >= 3){}
		else continue;
		if(isDebugCall) cout<<" Pass pt Cut ";

		// for counting to ensure things are not double counted
		if(passKin == false)
		{
			passKin = true;
			cutEventPass[cutFlowCH::Kinematics]++;
			cutPassHistW->Fill(cutFlowCH::Kinematics, weight);		
		}
		// Implement TriggerMatching
		InitTriggerMatchingTool();
		Bool_t passCutTrig = TriggerMatch( higgs_i->getLepton(), type);
		// For Trig eff studies
		if(anaType == doAnalysis::trigeff4l) passCutTrig = true;

		if(!passCutTrig) continue;
		if(isDebugCall) cout<<" Pass Trig ";
		
		// for counting to ensure things are not double counted
		 if(passTrig == false)
		{
			passTrig = true;
			cutEventPass[cutFlowCH::TriggerMatch]++;
			cutPassHistW->Fill(cutFlowCH::TriggerMatch, weight);					
		}
	
		// Find the closest one z mass
		Double_t tempDiff1 = fabs(higgs_i->getZ1()->get4Momentum()->M() - pdgZMass);
		Double_t tempDiff2 = fabs(higgs_i->getZ2()->get4Momentum()->M() - pdgZMass);

		// Debugging information
		if(isDebugCall)
		{
			cout<<" dM1: "<<setw(5)<< tempDiff1/1000;
			cout<<" dM2:"<<setw(5)<<tempDiff2/1000<<"  :Trigger Matched: "<<passCutTrig<<endl;
		}

		if(tempDiff1 < diffZ1Mass)
		{
			diffZ1Mass = tempDiff1;
			diffZ2Mass = tempDiff2; // added June 7, 2013 to fix the issue of chosing the right quadruplet
			temp = higgs_i;
		}
		else if(tempDiff1 == diffZ1Mass && tempDiff2 < diffZ2Mass)
		{
			diffZ2Mass = tempDiff2;
			temp = higgs_i;
		
		}
	}
	// Debugging information
	if(isDebugCall && temp != 0)
	{
		cout<<"Choose Quad Lepton pair: dM1: "<<setw(5)<< diffZ1Mass/1000<<" dM2:"<<setw(5)<<diffZ2Mass/1000<<endl;
		cout<<"Choose Quad Lepton pair: m4l: "<<setw(5)<< temp->get4Momentum()->M()/1000;
		cout<<" m12:"<<setw(5)<<temp->getZ1()->get4Momentum()->M()/1000;
		cout<<" m34:"<<setw(5)<<temp->getZ2()->get4Momentum()->M()/1000<<endl;
		cout<<"lepton chosen:"<<endl;
		Int_t i = 0;
		vector<ChargedLepton *> lep_i = temp->getLepton();
		for(vector<ChargedLepton *>::iterator itr_lep = lep_i.begin();
			itr_lep != lep_i.end(); ++itr_lep)
		{
			ChargedLepton* curr_i = *itr_lep;
			cout<<"Lepton: "<<i<<" pT:"<<curr_i->get4Momentum()->Pt()<<endl;
			i++;
		}
		lep_i.clear();
		cout<<"-------------------"<<endl;
	}
	else if(isDebugCall && temp == 0)
	{
		cout<<"Choose Quad Lepton pair: not intialized"<<endl;
	}


	return temp;
}
// Performs cuts on the quadlepton
Bool_t HiggsAnalysis::CutQuadLepton(QuadLepton * higgs, Int_t * cutEventPass, Int_t type, TH1D *cutPassHistW, Double_t *weightInit, Int_t * prodCH, Bool_t dolleeFlow)
{
	Bool_t passCut = true;
	// Getting the weight for each event
	if(doWeight && isMC) *weightInit = *weightInit * getHiggsWeight(higgs);
	Double_t weight = *weightInit;
	// Setting the weight for the plot histograms
	Hist->weight = weight * getCrossSectionWeight();

	// Z1 Mass cut
	Double_t candZ1Mass = higgs->getZ1()->get4Momentum()->M();
	Hist->Z1MassHist[type]->Fill(candZ1Mass, Hist->weight);
	if(candZ1Mass > 50*1000 && candZ1Mass < 106*1000) passCut = true;
	else passCut = false;
	// For Z->4l
	Bool_t passCutZ4l = true;
	if(candZ1Mass > 20*1000) passCutZ4l = true;
	else passCutZ4l = false;
	if(anaType == doAnalysis::Zllll) passCut = passCutZ4l;	
	if(!passCut) return passCut;
	cutEventPass[cutFlowCH::Z1Mass]++;
	cutPassHistW->Fill(cutFlowCH::Z1Mass, weight);

	// Z2 Mass Cut
	passCut = InterpolationZ2Cut(higgs->getZ2()->get4Momentum()->M(), higgs->get4Momentum()->M(), type);
	if(!passCut) return passCut;
	cutEventPass[cutFlowCH::Z2Mass]++;
	cutPassHistW->Fill(cutFlowCH::Z2Mass, weight);

	// Delta R Cut
	passCut = DeltaRCut(higgs->getLepton());
	if(!passCut) return passCut;
	
	// JPsi Veto
	passCut = JPsiVeto(higgs, type);
	if(!passCut) return passCut;		
	cutEventPass[cutFlowCH::DeltaR]++;
	cutPassHistW->Fill(cutFlowCH::DeltaR, weight);

	// Computing all the iso and d0sig
	vector<Double_t> pTCone20 = getTrackPT(higgs->getLepton());
	Bool_t passCutTrackIso = CutTrackIso(higgs->getLepton(), pTCone20, type, higgs);
	pTCone20.clear();

	vector<Double_t> eTCone20 = getCaloET(higgs->getLepton());
	Bool_t passCutCaloIso = CutCaloIso(higgs->getLepton(), eTCone20, type, higgs);
	eTCone20.clear();

	Bool_t passCutD0Sig = D0SigCut(higgs->getLepton(), type, higgs);

	// update the track iso if it passes track iso
	if(passCutTrackIso)
	{
		cutEventPass[cutFlowCH::TrackIso]++;
		cutPassHistW->Fill(cutFlowCH::TrackIso, weight);
	}
	// update the calo iso if it passes track iso & calo sip
	if(passCutTrackIso & passCutCaloIso)
	{
		cutEventPass[cutFlowCH::CaloIso]++;
		cutPassHistW->Fill(cutFlowCH::CaloIso, weight);
	}
	// update the track iso if it passes track iso
	if(passCutTrackIso & passCutCaloIso & passCutD0Sig)
	{
		cutEventPass[cutFlowCH::D0Sig]++;
		cutPassHistW->Fill(cutFlowCH::D0Sig, weight);
	}

	// For trigger effStudy
	if(anaType == doAnalysis::trigeff4l)
	{
		//Bool_t passCutTrig = false;
		//// For Counting DiMuon Trigger
		//if(type == quadType::Mu4) passCutTrig = SingleMuonTrigger(runNumber_sf) | DiMuonTrigger() ;

		//// for Counting electron trigger
		//if(type == quadType::El4) passCutTrig = SingleElectronTrigger(runNumber_sf)| DiElectronTrigger(runNumber_sf);

		//// If it passes any of the trigger, move ahead
		//if(type == quadType::Mu2El2 || type == quadType::El2Mu2)
		//{
		//	passCutTrig = SingleElectronTrigger(runNumber_sf)| DiElectronTrigger(runNumber_sf) | 
		//		SingleMuonTrigger(runNumber_sf) | DiMuonTrigger() |
		//		ElectronMuonTrigger();
		//}

		//if(!passCutTrig) return passCutTrig;

		//// Trigger match
		//passCutTrig = TriggerMatch( higgs->getLepton(), type);
		//if(!passCutTrig) return passCutTrig;
	}

	// Fill the output
	// Saving the weight
	higgs->weight_corr = weight;

	fillEventVarInfo(higgs, higgs->getQuadType(), prodCH);

	// if we are doing the llee cutflow, no need to store info for other CR as we might be overwriting. 
	if(dolleeFlow) return true;

	// Store the info for the CR before the iso and d0 cut
	if(doCR) outputTreeCR->fillTree(event, higgs,  higgs->getQuadType(), isMC);	

	// Final cuts
	passCut = passCutTrackIso & passCutCaloIso & passCutD0Sig;
	if(!passCut) return passCut;
 
	prodCH[higgs->getProductionChannel()]++;	
	// Fill the event var if it passes all cuts
	outputTree->fillTree(event, higgs,  higgs->getQuadType(), isMC);	
	if(doSysTree) outputTreeSys->fillTree(event, higgs,  higgs->getQuadType(), isMC);



	// Histrogram
	Double_t pTToStore [4] = {0};
	pTToStore[0] = higgs->getZ1()->getLepPlus()->get4Momentum()->Pt();
	pTToStore[1] = higgs->getZ1()->getLepNeg()->get4Momentum()->Pt();
	pTToStore[2] = higgs->getZ2()->getLepPlus()->get4Momentum()->Pt();
	pTToStore[3] = higgs->getZ2()->getLepNeg()->get4Momentum()->Pt();	
	std::sort(pTToStore, pTToStore + 4);
	Hist->quadFinalPtHist[type][0]->Fill(pTToStore[0], Hist->weight);
	Hist->quadFinalPtHist[type][1]->Fill(pTToStore[1], Hist->weight);
	Hist->quadFinalPtHist[type][2]->Fill(pTToStore[2], Hist->weight);
	Hist->quadFinalPtHist[type][3]->Fill(pTToStore[3], Hist->weight);

	// For Debugging 
	if(isDebugCall)
	{
		cout<<"----------------------"<<endl;
		cout<<"Leptons in the Quadruplet"<<endl;

		vector <ChargedLepton *> lepList = higgs->getLepton();
		Int_t i = 0;
		for(vector<ChargedLepton *>::iterator itr_i = lepList.begin();
		itr_i != lepList.end(); ++itr_i)
		{
			// Getting the charged Lepton
			ChargedLepton * lep_i = *itr_i;
			if (lep_i->getFlavor() == flavor::Electron)
			{
				cout<<"Lepton "<< i <<" Type:Electron ";
				cout<<"Eta: "<<lep_i->get4Momentum()->Eta();
				cout<<" Pt:"<<lep_i->get4Momentum()->Pt()<<endl;
			}
			if (lep_i->getFlavor() == flavor::Muon)
			{
				cout<<"Lepton "<< i <<" Type:Muon";
				if(lep_i->getType() == leptonType::MuonCalo) cout<<"calo ";
				else if(lep_i->getType() == leptonType::MuonStaco) cout<<"Staco ";
				else if(lep_i->getType() == leptonType::MuonStandAlone) cout<<"StandAlone ";
				
				cout<<"Eta: "<<lep_i->get4Momentum()->Eta();
				cout<<" Pt:"<<lep_i->get4Momentum()->Pt()<<endl;
			}
			i++;
		}
		lepList.clear();
		cout<<"----------------------"<<endl;
		cout<<"----------------------"<<endl<<endl<<endl;
	}
	
	return true;
}
// Helper for Z2 mass cut
Bool_t HiggsAnalysis::InterpolationZ2Cut(Double_t Z2Mass, Double_t QuadMass, Int_t type)
{
	Bool_t passCut = false;

	const Int_t nBin = 2;
	Double_t mass[nBin] = {140.0*1000, 190.0*1000};
	Double_t cut[nBin] = {12.0*1000, 50.0*1000};

	// For Z->4l Analysis overwrite the cuts above
	if(anaType == doAnalysis::Zllll)
	{
		mass[0] = -1;
		mass[1] = 999999999;

		cut[0] = 5.0 * 1000;
		cut[1] = 5.0 * 1000;
	}
	
	Double_t z2CutVal = 0;
	Int_t index = -1;

	for(Int_t i = 0; i < nBin; i++)
	{
		if(QuadMass > mass[i]) index = i;
	}

	if(index == -1) z2CutVal = cut[0];
	else if (index == (nBin - 1)) z2CutVal = cut[(nBin - 1)];
	else z2CutVal = cut[index] + (QuadMass - mass[index]) * (cut[index+1] - cut[index])/(mass[index+1] - mass[index]);

	Hist->Z2M4lHist[type]->Fill(Z2Mass, QuadMass, Hist->weight);

	if(Z2Mass > z2CutVal && Z2Mass < 115*1000) passCut = true;
	
	return passCut;
}
// Delta R cuts
Bool_t HiggsAnalysis::DeltaRCut(vector<ChargedLepton *> curr_lep)
{
	Bool_t passCut = true;
	for(vector<ChargedLepton *>::iterator itr_i = curr_lep.begin();
		itr_i != curr_lep.end(); ++itr_i)
	{
		for(vector<ChargedLepton *>::iterator itr_j = itr_i+1;
			itr_j != curr_lep.end(); ++itr_j)
		{
			// Getting the charged Lepton
			ChargedLepton * lep_i = *itr_i;
			ChargedLepton * lep_j = *itr_j;
		
			// Since the cut value is dependant on flavour type
			// getting that
			Double_t cutDeltaR = 0;
			if(lep_i->getFlavor() == lep_j->getFlavor()) {cutDeltaR = 0.10;}
			else {cutDeltaR = 0.20;}

			// Actually compute the delta R and cut on it
			Double_t dataDeltaR = 0;
			dataDeltaR = lep_i->get4Momentum()->DeltaR(*(lep_j->get4Momentum()));
			
			// Debugging information
			if(isDebugCall){cout<<"Delta R "<<dataDeltaR<<"\t";}

			if(dataDeltaR < cutDeltaR) passCut = false;

			if(!passCut) break;
		}
			// Debugging information
			if(isDebugCall){cout<<endl;}

		if(!passCut) break;
	}
	return passCut;
}

// JPsi Veto
Bool_t HiggsAnalysis::JPsiVeto(QuadLepton * higgs, Int_t type)
{
	Bool_t passCut = true;
	if(type == analysisType::El2Mu2 ||type == analysisType::Mu2El2) return true;

	TLorentzVector crossPair1 = *(higgs->getZ1()->getLepPlus()->get4Momentum()) +  
								*(higgs->getZ2()->getLepNeg()->get4Momentum());
	TLorentzVector crossPair2 = *(higgs->getZ2()->getLepPlus()->get4Momentum()) +  
								*(higgs->getZ1()->getLepNeg()->get4Momentum());
	// Histrogram
	Hist->JPsiHist[type]->Fill(crossPair1.M(), Hist->weight);
	Hist->JPsiHist[type]->Fill(crossPair2.M(), Hist->weight);
	
	if(crossPair1.M() < 5 * 1000 || crossPair2.M() < 5* 1000) { passCut = false;}

	return passCut;
}

// Getting the pT weights for each lepton
vector<Double_t> HiggsAnalysis::getTrackPT(vector<ChargedLepton *> curr_lep)
{
	// Container
	vector <Double_t> trackPt;

	// Delta R cut
	Double_t cutDeltaR = 0.20;

	vector<ChargedLepton *>::iterator itr_i;
	vector<ChargedLepton *>::iterator itr_j;
	int i = 0;
	for(itr_i = curr_lep.begin(); itr_i != curr_lep.end(); itr_i++)
	{
		i++;
		ChargedLepton *lep_i = *itr_i;
		Double_t curr_PTCone20 = lep_i->getPTCone20();
		if(isDebugCall) cout<<"-------------------"<<endl<<"Lepton Track Iso Original PTCone20: "<<curr_PTCone20<<endl;
		
		for(itr_j = curr_lep.begin(); itr_j != curr_lep.end(); itr_j++)
		{
			ChargedLepton *lep_j = *itr_j;
			
			// Don't compare the same thing
			if(lep_i == lep_j) continue;

			Double_t dataDeltaR = lep_i->get4Momentum()->DeltaR(*(lep_j->get4Momentum()));
			Double_t toSubtract = 0;
			if(isDebugCall) cout<<"lepton "<< i <<" delta R "<<dataDeltaR<<endl;			
			if(dataDeltaR < cutDeltaR) 
			{
				if(lep_j->getType() == leptonType::MuonStaco || lep_j->getType() == leptonType::MuonCalo)
				{
					toSubtract = (1/fabs((lep_j->GetMuon()->id_qoverp()))*sin((lep_j->GetMuon()->id_theta())));
					if(isDebugCall) cout<<"Muon "<< i <<" ToSubtract "<<toSubtract<<endl;
					
				}
				else if (lep_j->getType() == leptonType::MuonStandAlone)
				{
					toSubtract = 0;
				}
				else if(lep_j->getType() == leptonType::ElectronGSF)
				{
					toSubtract = lep_j->GetElectron()->trackpt();
					if(isDebugCall) cout<<"Electron "<< i <<" ToSubtract "<<toSubtract<<endl;
				}
				else
				{
					cout<<"UnExpected lepton:: trackIso Cut"<<endl;
				}
			}
			curr_PTCone20 = curr_PTCone20 - toSubtract;
		}
		if(isDebugCall) cout<<"Track Iso Final PTCone20: "<<curr_PTCone20<<endl;
		trackPt.push_back(curr_PTCone20);
	}

	return trackPt;
	
}
// Does the track iso cut
Bool_t HiggsAnalysis::CutTrackIso(vector<ChargedLepton *> curr_lep, vector<Double_t> trackPT, Int_t type, QuadLepton* higgs)
{
 
	higgs->trackIso.clear();

	Double_t cutTrackIso = 0.15;

	Bool_t cutAll = true;
	for(UInt_t i = 0; i < curr_lep.size(); i++)
	{
		ChargedLepton * lep_i = curr_lep[i];
		Double_t dataTrackIso = 0;
		dataTrackIso = trackPT[i]/lep_i->get4Momentum()->Pt();
		if(isDebugCall) cout<<"Track Iso: "<<dataTrackIso<<endl;
		// histrogram
		Hist->TrackIsoHist[type]->Fill(dataTrackIso, Hist->weight);

		// Store track iso
		higgs->valTrackIso.push_back(dataTrackIso);

		if(dataTrackIso > cutTrackIso) higgs->trackIso.push_back(false);
		else higgs->trackIso.push_back(true);	
	}

	// Sanity check
	if((Int_t) higgs->trackIso.size() != 4) cout<<"ERROR: Trackiso vector not size"<<endl;

	for(Int_t i = 0; i < (Int_t) higgs->trackIso.size(); i++) cutAll = cutAll & higgs->trackIso[i];

	return cutAll;
}
// Gets the calo isolation weight. Currently no calo isolation correction
vector<Double_t> HiggsAnalysis::getCaloET(vector<ChargedLepton *> curr_lep)
{
	// Container
	vector <Double_t> caloEt;

	//Double_t *EtOverlap = new Double_t[curr_lep.size()];
	//Double_t *EtCone20CorrectedIso =  new Double_t[curr_lep.size()];
	Double_t EtOverlap[4]; // = new Double_t[curr_lep.size()];
	Double_t EtCone20CorrectedIso[4];// =  new Double_t[curr_lep.size()];
	vector<ChargedLepton *>::iterator itr_i;
	vector<ChargedLepton *>::iterator itr_j;

	Double_t deltaRCut = 0.18;

	Int_t i = -1;
	// This gets what I need to subtract in the end for calo iso
	for(itr_i = curr_lep.begin(); itr_i != curr_lep.end(); itr_i++)
	{
		// Just to keep track of what event loop we are at
		i++;

		ChargedLepton *lep_i = *itr_i;

		if(lep_i->getFlavor() == flavor::Electron)
		{
			D3PDReader::ElectronD3PDObjectElement* el_curr = lep_i->GetElectron();
			//if(dataYear == 2011)
			//{
			//	EtOverlap[i] = el_curr->cl_E()/cosh(el_curr->tracketa());
			//}
			//else if(dataYear == 2012)
			//{
			//	EtOverlap[i] = el_curr->rawcl_E()/cosh(el_curr->tracketa());
			//}
			EtOverlap[i] = el_curr->rawcl_E()/cosh(el_curr->tracketa());
		}
		else if(lep_i->getFlavor() == flavor::Muon)
		{
			D3PDReader::MuonD3PDObjectElement* mu_curr = lep_i->GetMuon();
			EtOverlap[i] = mu_curr->pt();
			EtOverlap[i] = 0;
		}
	}
	// This gets what to actually subtract from
	// This is were the calo isolation corrections will be implemented
	i = -1;
	if(isDebugCall) cout<<"Calo EtCone20CorrectedIso: "<<endl;
	for(itr_i = curr_lep.begin(); itr_i != curr_lep.end(); itr_i++)
	{
		// Just to keep track of what event loop we are at
		i++;

		ChargedLepton *lep_i = *itr_i;

		if(lep_i->getFlavor() == flavor::Electron)
		{
			D3PDReader::ElectronD3PDObjectElement* el_curr = lep_i->GetElectron();
			if(dataYear == 2011)
			{
				EtCone20CorrectedIso[i] = el_curr->Etcone20();
			}
			else if(dataYear == 2012)
			{
				EtCone20CorrectedIso[i] = el_curr->topoEtcone20();
			}
			// Smearing Correction
			if(doCaloIsolationCorr && doCorr) EtCone20CorrectedIso[i] = getCaloIsoCorrection(el_curr);
		}
		else if(lep_i->getFlavor() == flavor::Muon)
		{
			D3PDReader::MuonD3PDObjectElement* mu_curr = lep_i->GetMuon();
			EtCone20CorrectedIso[i] = mu_curr->etcone20();
		}
		if(isDebugCall) cout<<"Lepton i: "<<i<<" "<< EtCone20CorrectedIso[i]<<endl;
		
	}
	// Now to actually do the subtraction
	i = -1; 
	if(isDebugCall) cout<<"Final after subraction: "<<endl;	
	for(itr_i = curr_lep.begin(); itr_i != curr_lep.end(); itr_i++)
	{
		// Just to keep track of what event loop we are at
		i++;
		ChargedLepton *lep_i = *itr_i;
		
		Int_t j = -1;
		for(itr_j = curr_lep.begin(); itr_j != curr_lep.end(); itr_j++)
		{
			// Just to keep track of what event loop we are at
			j++;
			ChargedLepton *lep_j = *itr_j;
			// Don't want to comapare to the same thing
			if(lep_i == lep_j) continue;

			if(DeltaRCaloIso(lep_i, lep_j) < deltaRCut)
			{
				EtCone20CorrectedIso[i] = EtCone20CorrectedIso[i] - EtOverlap[j];
			}
		}
		caloEt.push_back(EtCone20CorrectedIso[i]);
		if(isDebugCall) cout<<"Lepton i: "<<i<<" "<< EtCone20CorrectedIso[i]<<endl;
		
	}
	return caloEt;
}
// For CaloIso Corrections
Double_t HiggsAnalysis::getCaloIsoCorrection(D3PDReader::ElectronD3PDObjectElement* el_curr)
{
	Double_t caloIsoCorr = -1;
	if(dataYear == 2011 && doCaloIsolationCorr)
	{
		caloIsoCorr = CaloIsoCorrection::GetNPVCorrectedIsolation( 	getNVertex(2),
																	el_curr->etas2(),
																	20, 
  																	isMC,
																	el_curr->Etcone20(),
                                                             		CaloIsoCorrection::ELECTRON,
                                                             		CaloIsoCorrection::REL17);
	}
	else if(dataYear == 2012 && doCaloIsolationCorr)
	{
		Double_t clusterE = el_curr->cl_E();
		if(doCorr) clusterE = el_curr->cl_E_unsmeared;
		caloIsoCorr = CaloIsoCorrection::GetPtEDCorrectedTopoIsolation(	el_curr->ED_median(),
																		clusterE,
																		el_curr->etas2(),
																		el_curr->etap(),
																		el_curr->cl_eta(),
                                                                    	20,
                                                                    	isMC,
																		el_curr->topoEtcone20(),
                                                                    	false, 
                                                                    	CaloIsoCorrection::ELECTRON,
                                                                    	CaloIsoCorrection::REL17_2);
	}
	return caloIsoCorr;
}


// Cuts on Calo Isolaton
Bool_t HiggsAnalysis::CutCaloIso(vector<ChargedLepton *> curr_lep, vector<Double_t> caloPT, Int_t type, QuadLepton* higgs)
{
	vector<ChargedLepton *>::iterator itr_i;
	// Counter
	Int_t i = -1;
	higgs->caloIso.clear();
	
	// Cut Values;
	Double_t muonCut = 0.30;
	Double_t muonStandAloneCut = 0.15;

	Double_t electronCut2011 = 0.30;
	Double_t electronCut2012 = 0.20;

	Bool_t cutAll = true;

	// Actual Cuts
	for(itr_i = curr_lep.begin(); itr_i != curr_lep.end(); itr_i++)
	{
		i++;
		ChargedLepton *lep_i = *itr_i;

		Double_t dataCut = caloPT[i]/lep_i->get4Momentum()->Pt();
		
		// Store track iso
		higgs->valCaloIso.push_back(dataCut);	
		
		Bool_t cutLep = true;

		// Debugging information
		if(isDebugCall)
		{
			cout<<"-------------------"<<endl<<"Calo Iso: ";
			if(lep_i->getFlavor() == flavor::Electron) {cout<<"electron: ";}
			else if(lep_i->getType() == leptonType::MuonStandAlone) {cout<<"StandAloneMuon: ";}
			else if(lep_i->getType() == leptonType::MuonCalo) {cout<<"CaloMuon: ";}
			else if(lep_i->getType() == leptonType::MuonStaco) {cout<<"StacoMuon: ";}
			
			cout<<dataCut << endl;
		}

		if(lep_i->getFlavor() == flavor::Electron)
		{
			// histrogram
			Hist->elCaloIsoHist[type]->Fill(dataCut, Hist->weight);
			if(dataYear == 2011) 
			{
				if(dataCut > electronCut2011) cutLep = false;
			}
			else if(dataYear == 2012)
			{
				if(dataCut > electronCut2012) cutLep = false;
			}
		}
		if(lep_i->getFlavor() == flavor::Muon)
		{
			// histrogram
			Hist->muCaloIsoHist[type]->Fill(dataCut, Hist->weight);

			if(lep_i->getType() != leptonType::MuonStandAlone) 
			{
				if(dataCut > muonCut) cutLep = false;
			}
			else if(lep_i->getType() == leptonType::MuonStandAlone)
			{
				if(dataCut > muonStandAloneCut) cutLep = false;
			}
		}
		higgs->caloIso.push_back(cutLep);
	}
	// Sanity check
	if((Int_t) higgs->caloIso.size() != 4) cout<<"ERROR: caloiso vector not size"<<endl;

	for(Int_t i = 0; i < (Int_t) higgs->caloIso.size(); i++) cutAll = cutAll & higgs->caloIso[i];

	return cutAll;
}

Bool_t HiggsAnalysis::D0SigCut(vector<ChargedLepton *> curr_lep, Int_t type, QuadLepton* higgs)
{
	vector<ChargedLepton *>::iterator itr_i;
	higgs->d0Sig.clear();
	// Cut Values;
	Double_t muD0SigCut = 3.5;
	Double_t elD0SigCut = 6.5;

	Bool_t cutAll = true;

	if(isDebugCall) cout<<"-------------------"<<endl<<"D0 Sig:"<<endl;
	// Actual Cuts
	for(itr_i = curr_lep.begin(); itr_i != curr_lep.end(); itr_i++)
	{
		ChargedLepton *lep_i = *itr_i;

		Bool_t cutLep = true;

		if(lep_i->getFlavor() == flavor::Electron)
		{
			Double_t dataD0Sig = fabs(lep_i->GetElectron()->trackd0pvunbiased())/ lep_i->GetElectron()->tracksigd0pvunbiased();

			higgs->valD0Sig.push_back(dataD0Sig);	
			
			// histrogram
			Hist->elSigHist[type]->Fill(dataD0Sig, Hist->weight);
			if(isDebugCall) cout<<"Electron d0Sig: "<<dataD0Sig<<endl;				
			if(dataD0Sig > elD0SigCut) cutLep = false;
		}
		if(lep_i->getFlavor() == flavor::Muon)
		{
			Double_t dataD0Sig = fabs(lep_i->GetMuon()->trackd0pvunbiased())/ lep_i->GetMuon()->tracksigd0pvunbiased();
			
			higgs->valD0Sig.push_back(dataD0Sig);	
			
			// histrogram
			Hist->muSigHist[type]->Fill(dataD0Sig, Hist->weight);
			if(isDebugCall) cout<<"Muon d0Sig: "<<dataD0Sig<<endl;	
			if(dataD0Sig > muD0SigCut) cutLep = false;			
		}
		higgs->d0Sig.push_back(cutLep);
	}
	// Sanity check
	if((Int_t) higgs->d0Sig.size() != 4) cout<<"ERROR: d0Sig vector not size"<<endl;

	for(Int_t i = 0; i < (Int_t) higgs->d0Sig.size(); i++) cutAll = cutAll & higgs->d0Sig[i];

	return cutAll;


}

// Trigger Matching
Bool_t HiggsAnalysis::TriggerMatch(vector<ChargedLepton *> curr_lep ,Int_t type)
{

	// Getting the Trigger Strings
	TString singleEl [2];
	TString diEl [2];
	TString singleMu [2];
	TString diMu [2];
	TString eMu [2];
 	FillTriggerString(singleMu, diMu, singleEl, diEl, eMu);

	// Boolean var to keep track
	Bool_t muTrigger = false;
	Bool_t elTrigger = false;
	Bool_t emuTrigger = false;

	if(type == analysisType::Mu4 || type == analysisType::El2Mu2 || type == analysisType::Mu2El2)
	{
		Bool_t singleMuTrigger = TriggerMatchSingleMuon(curr_lep, singleMu);
		Bool_t diMuTrigger = TriggerMatchDiMuon(curr_lep, diMu);
		muTrigger = (singleMuTrigger || diMuTrigger);
	}
	if(type == analysisType::El4 || type == analysisType::El2Mu2 || type == analysisType::Mu2El2)
	{
		Bool_t singleElTrigger = TriggerMatchSingleElectron(curr_lep, singleEl); 
		Bool_t diElTrigger = TriggerMatchDiElectron(curr_lep, diEl); 
		elTrigger = (singleElTrigger || diElTrigger);
	}
	if(type == analysisType::El2Mu2 || type == analysisType::Mu2El2)
	{
		Bool_t intEMuTrigger = TriggerMatchElectronMuon(curr_lep, eMu); 
		emuTrigger = intEMuTrigger;
	}

	if(type == analysisType::Mu4) return muTrigger;
	if (type == analysisType::El4) return elTrigger;
	if (type == analysisType::El2Mu2 || type == analysisType::Mu2El2) return (muTrigger || elTrigger || emuTrigger);
	else{ cout <<"Error for type in trigger matching"<<endl;}
	
	return true;	
}
// Single Muon
Bool_t HiggsAnalysis::TriggerMatchSingleMuon(vector<ChargedLepton *> curr_lep, TString triggerName[])
{
	Bool_t passCutTrig = false;
	Bool_t passCutThres = false;

	// Threshold
	Double_t pTThreshold = -1*1000;

	// Some random fix... got it from Fabien
	if (D3DPTriggerDev){
 		triggerNavigationVariables->set_trig_RoI_EF_mu_TrigMuonEFInfoContainerStatus(event->trig_RoI_EF_mu.TrigMuonEFInfoContainer_eMuonEFInfoStatus());
        triggerNavigationVariables->set_trig_RoI_EF_mu_TrigMuonEFInfoContainer(event->trig_RoI_EF_mu.TrigMuonEFInfoContainer_eMuonEFInfo());
	}  

	// Looking to see it the SingleMuon Trigger was flagged or not
	Bool_t singleMuonTrigged = SingleMuonTrigger(runNumber_sf);
	// Return if the singleMuon Trigger was not flagged
	if(!singleMuonTrigged) return false;

	// Now check if any of the leptons have triggered the SingleMuon trigger
	// Iterators
	vector<ChargedLepton *>::iterator itr_i;

	for(itr_i = curr_lep.begin(); itr_i != curr_lep.end(); itr_i++)
	{
		ChargedLepton *lep_i = *itr_i;
		// Don't want to compare an electron...
		if(lep_i->getFlavor() != flavor::Muon) continue;

		// No Trigger for Calo Muon
		if(lep_i->getType() == leptonType::MuonCalo) continue;

		for(Int_t i = 0; i < 2; i++)
		{
			// If the second trigger is not initialized
			if(triggerName[i] == "") continue;

			Bool_t trigMatch = muonTriggerMatchTool->match(lep_i->GetMuon()->eta(), 
            	                            		 	   lep_i->GetMuon()->phi(), 
                	                         			 	triggerName[i].Data());

			if(trigMatch == true) passCutTrig = true;

			// For Threshold Cuts currently non-existance
			if(passCutTrig)
			{
				if(lep_i->GetMuon()->pt() > pTThreshold){passCutThres = true;}
			}
			
		}

	}
	
	return passCutThres;
}
// Di Muons Trigger matching
Bool_t HiggsAnalysis::TriggerMatchDiMuon(vector<ChargedLepton *> curr_lep, TString triggerName[])
{
	Bool_t passCutTrig = false;
	Bool_t passCutThres = false;

	// Threshold
	Double_t pTThreshold1 = -1*1000;
	Double_t pTThreshold2 = -1*1000;

	// Some random fix... got it from Fabien
	if (D3DPTriggerDev){
 		triggerNavigationVariables->set_trig_RoI_EF_mu_TrigMuonEFInfoContainerStatus(event->trig_RoI_EF_mu.TrigMuonEFInfoContainer_MuonEFInfoStatus());
        triggerNavigationVariables->set_trig_RoI_EF_mu_TrigMuonEFInfoContainer(event->trig_RoI_EF_mu.TrigMuonEFInfoContainer_MuonEFInfo());
	}  

	// Looking to see it the SingleMuon Trigger was flagged or not
	Bool_t diMuonTrigged = DiMuonTrigger();
	// Return if the diMuon Trigger was not flagged
	if(!diMuonTrigged) return false;

	// Now check if any of the leptons have triggered the SingleMuon trigger
	// Iterators
	vector<ChargedLepton *>::iterator itr_i;
	vector<ChargedLepton *>::iterator itr_j;

	// To store the result from the matching tool
	std::pair<bool, bool> PassT1, PassT2;
	std::pair<bool, bool> PassT3, PassT4;
	// Looping over all the events
	for(itr_i = curr_lep.begin(); itr_i != curr_lep.end(); itr_i++)
	{
		ChargedLepton *lep_i = *itr_i;
		// Don't want to compare an electron...
		if(lep_i->getFlavor() != flavor::Muon) continue;
		// No Trigger for Calo Muon
		if(lep_i->getType() == leptonType::MuonCalo) continue;

		for(itr_j = curr_lep.begin(); itr_j != curr_lep.end(); itr_j++)
		{
			ChargedLepton *lep_j = *itr_j;
			// Don't want to compare to the same thing
			if(lep_i == lep_j) continue;
			// Don't want to compare an electron...
			if(lep_j->getFlavor() != flavor::Muon) continue;
			// No Trigger for Calo Muon
			if(lep_j->getType() == leptonType::MuonCalo) continue;

			for(Int_t i = 0; i < 2; i++)
			{
				// If the second trigger is not initialized
				if(triggerName[i] == "") continue;

				Bool_t trigMatch1 = muonTriggerMatchTool->matchDimuon(*(lep_i->get4Momentum()), 
            		                            		 	  		 *(lep_j->get4Momentum()), 
                		                         			 		  triggerName[i].Data(), 
																	  PassT1, PassT2);
				Bool_t trigMatch2 = muonTriggerMatchTool->matchDimuon(*(lep_j->get4Momentum()), 
            		                            		 	  		 *(lep_i->get4Momentum()), 
                		                         			 		  triggerName[i].Data(), 
																	  PassT3, PassT4);

			
				if(!(trigMatch1||trigMatch2)) cout<<"DiMuon Error: Trigger not supported: "<<triggerName[i]<<endl;
				else
				{
        			passCutTrig=(PassT1.first && PassT2.second)||(PassT3.first && PassT4.second);
				}
				// For Threshold Cuts currently non-existance
				if(passCutTrig)
				{
					if(lep_i->GetMuon()->pt() > pTThreshold1 && lep_j->GetMuon()->pt() > pTThreshold2 ){passCutThres = true;}
				}
			}
		}
	}
	return passCutThres;
}
// Single Electron
Bool_t HiggsAnalysis::TriggerMatchSingleElectron(vector<ChargedLepton *> curr_lep, TString triggerName[])
{
	Bool_t passCutTrig = false;
	Bool_t passCutThres = false;

	// Threshold
	Double_t pTThreshold = -1*1000;

	// Looking to see it the SingleMuon Trigger was flagged or not
	Bool_t singleElecronTrigged = SingleElectronTrigger(runNumber_sf);
	// Return if the singleMuon Trigger was not flagged
	if(!singleElecronTrigged) return false;

	// Now check if any of the leptons have triggered the SingleMuon trigger
	// Iterators
	vector<ChargedLepton *>::iterator itr_i;

	for(itr_i = curr_lep.begin(); itr_i != curr_lep.end(); itr_i++)
	{
		ChargedLepton *lep_i = *itr_i;
		// Don't want to compare an electron...
		if(lep_i->getFlavor() != flavor::Electron) continue;

		for(Int_t i = 0; i < 2; i++)
		{
			// If the second trigger is not initialized
			if(triggerName[i] == "") continue;

			Bool_t trigMatch = 
				electronTriggerMatchTool->match(lep_i->GetElectron()->tracketa(), 
                                          		lep_i->GetElectron()->trackphi(), 
                                          		 triggerName[i].Data());

			
			if(trigMatch == true) passCutTrig = true;
			
			// For Threshold Cuts currently non-existant for 2012 Moriond	
			if(passCutTrig)
			{ 
				if(lep_i->GetElectron()->pt() > pTThreshold){passCutThres = true;}
			}
		}
	}
	
	return passCutThres;
}
// Di Electron trigger matching
Bool_t HiggsAnalysis::TriggerMatchDiElectron(vector<ChargedLepton *> curr_lep, TString triggerName[])
{
	Bool_t passCutTrig = false;
	Bool_t passCutThres = false;

	// Threshold
	Double_t pTThreshold1 = -1*1000;
	Double_t pTThreshold2 = -1*1000;

	// Looking to see it the diElectron Trigger was flagged or not
	Bool_t diElectronTrigged = DiElectronTrigger(runNumber_sf);

	// Return if the diMuon Trigger was not flagged
	if(!diElectronTrigged) return false;

	// Now check if any of the leptons have triggered the SingleMuon trigger
	// Iterators
	vector<ChargedLepton *>::iterator itr_i;
	vector<ChargedLepton *>::iterator itr_j;

	// To store the result from the matching tool
	Bool_t PassT1, PassT2;

	// Looping over all the events
	for(itr_i = curr_lep.begin(); itr_i != curr_lep.end(); itr_i++)
	{
		ChargedLepton *lep_i = *itr_i;
		// Don't want to compare an muon...
		if(lep_i->getFlavor() != flavor::Electron) continue;

		for(itr_j = itr_i+1; itr_j != curr_lep.end(); itr_j++)
		{
			ChargedLepton *lep_j = *itr_j;
			// Don't want to compare to the same thing
			if(lep_i == lep_j) continue;
			// Don't want to compare an electron...
			if(lep_j->getFlavor() != flavor::Electron) continue;
	
			for(Int_t i = 0; i < 2; i++)
			{
				// If the second trigger is not initialized
				if(triggerName[i] == "") continue;

				Bool_t trigMatch = electronTriggerMatchTool->matchDielectron(*(lep_i->get4Momentum()), 
            		                            		 	  		 *(lep_j->get4Momentum()), 
                		                         			 		  triggerName[i].Data(), 
																	  PassT1, PassT2);
		
				if(!trigMatch) cout<<"DiElectron Error: Trigger not supported: "<<triggerName[i]<<endl;
				else
				{
        			passCutTrig=(PassT1 && PassT2);
				}
				// For Threshold Cuts currently non-existant for 2012 Moriond				
				if(passCutTrig)
				{
					if(lep_i->GetElectron()->pt() > pTThreshold1 && lep_j->GetElectron()->pt() > pTThreshold2 ){ passCutThres = true; }
				}
			}
		}
	}

	return passCutThres;
}
// Emu trigger matching
Bool_t HiggsAnalysis::TriggerMatchElectronMuon(vector<ChargedLepton *> curr_lep, TString triggerName[])
{
	Bool_t passCutTrig = false;
	Bool_t passCutThres = false;

	// Threshold
	Double_t pTThreshold1 = -1*1000;
	Double_t pTThreshold2 = -1*1000;

	// Some random fix... got it from Fabien
	if (D3DPTriggerDev){
 		triggerNavigationVariables->set_trig_RoI_EF_mu_TrigMuonEFInfoContainerStatus(event->trig_RoI_EF_mu.TrigMuonEFInfoContainer_MuonEFInfoStatus());
        triggerNavigationVariables->set_trig_RoI_EF_mu_TrigMuonEFInfoContainer(event->trig_RoI_EF_mu.TrigMuonEFInfoContainer_MuonEFInfo());
	}  

	// Looking to see it the SingleMuon Trigger was flagged or not
	Bool_t electronMuonTrigged = ElectronMuonTrigger();
	// Return if the diMuon Trigger was not flagged
	if(!electronMuonTrigged) return false;

	// Now check if any of the leptons have triggered the SingleMuon trigger
	// Iterators
	vector<ChargedLepton *>::iterator itr_i;
	vector<ChargedLepton *>::iterator itr_j;

	// Looping over all the events
	for(itr_i = curr_lep.begin(); itr_i != curr_lep.end(); itr_i++)
	{
		ChargedLepton *lep_i = *itr_i;
		// First get the electron
		if(lep_i->getFlavor() != flavor::Electron) continue;

		for(itr_j = curr_lep.begin(); itr_j != curr_lep.end(); itr_j++)
		{
			ChargedLepton *lep_j = *itr_j;
			// Don't want to compare to the same thing
			if(lep_i == lep_j) continue;
			// Get the muon
			if(lep_j->getFlavor() != flavor::Muon) continue;
	
			for(Int_t i = 0; i < 2; i++)
			{
				// If the second trigger is not initialized
				if(triggerName[i] == "") continue;

				Bool_t trigMatch = electronTriggerMatchTool->matchElectronMuon(*(lep_i->get4Momentum()), 
            		                            		 	  		 *(lep_j->get4Momentum()), 
                		                         			 		  triggerName[i].Data());

				if(trigMatch) 
				{
        			passCutTrig=true;;
				}
				// For Threshold Cuts currently non-existant for 2012 Moriond					
				if(passCutTrig)
				{
					if(lep_i->GetElectron()->pt() > pTThreshold1 && lep_j->GetMuon()->pt() > pTThreshold2 )
					{ passCutThres = true;}
				}
			}
		}
	}
	return passCutThres;
}
// Find the mass and its uncertainity...
// Also Does FSR and Z mass Constraint
void HiggsAnalysis::MassCalc(QuadLepton * higgs, Int_t Type)
{
	// Setting the mass
	Double_t mass = InvariantMass::invMassCalc(higgs->getLeptonLorentz());
	higgs->SetMass(mass);
	higgs->SetElRescale(corr->elRescale);
	higgs->FillCovMatrix(runNumber_sf);
	// Getting the Mass Error
	Double_t massErr = -1;
	massErr = MassError::sigmaMassCalc_d0z0PhiThetaP(higgs->getLeptonLorentz(),higgs->getLeptonCovMatrix());
	higgs->SetMassErr(massErr);

	higgs->massID = InvariantMass::invMassCalc(higgs->getLeptonLorentz(muonType::ID));
	higgs->massMS = InvariantMass::invMassCalc(higgs->getLeptonLorentz(muonType::MS));
	higgs->massErrID = MassError::sigmaMassCalc_d0z0PhiThetaP(higgs->getLeptonLorentz(muonType::ID),higgs->getLeptonCovMatrix(muonType::ID));
	higgs->massErrMS = MassError::sigmaMassCalc_d0z0PhiThetaP(higgs->getLeptonLorentz(muonType::MS),higgs->getLeptonCovMatrix(muonType::MS));
	
	// Creating a sum of TLorentzVector for the overall sum
	vector<TLorentzVector> leptonLorentz = higgs->getLeptonLorentz();
	for(Int_t i = 0; i < (Int_t) leptonLorentz.size(); i++)
	{
		higgs->sum_unconstrained = higgs->sum_unconstrained + leptonLorentz[i];
	}
	leptonLorentz.clear();

	// Actual FSR
	CorrectFSR(higgs, Type);
	
	// Z Mass Constraint Fix
	CorrectZMassConstraint(higgs);
	CorrectZMassConstraint(higgs, muonType::ID);
	CorrectZMassConstraint(higgs, muonType::MS);
	return;
}
// FSR Corrections
void HiggsAnalysis::CorrectFSR(QuadLepton * higgs, Int_t Type)
{
	if (curEvent == 6 || curEvent == 73) isDebugCall = true;
	// Default FSR correction - no FSR correction
	higgs->fsrType = fsrType::noFSR;
	// Getting all the lepton in the quad
	vector<ChargedLepton *> lep_cont = higgs->getLepton();

	vector<FsrPhotons::FsrCandidate> candVec;

	// Calculating the iso for the photon
	//
	vector<float>* ph_Etcone40_corr = new vector<float>();
	double caloIsoCorr = 0;
	for(Int_t i = 0; i < event->ph.n(); i++)
	{
		if(dataYear == 2011 )
		{
			caloIsoCorr = CaloIsoCorrection::GetNPVCorrectedIsolation( 	getNVertex(2),
																		event->ph[i].etas2(),
																		40, 
  																		isMC,
																		event->ph[i].Etcone40(),
    	                                                         		CaloIsoCorrection::PHOTON,
    	                                                         		CaloIsoCorrection::REL17);
			ph_Etcone40_corr->push_back(caloIsoCorr);
		}
		else if(dataYear == 2012)
		{
			Double_t clusterE = event->ph[i].cl_E();
			if(doCorr) clusterE = event->ph[i].cl_E_unsmeared;
			caloIsoCorr = CaloIsoCorrection::GetPtEDCorrectedTopoIsolation(	event->ph[i].ED_median(),
																			clusterE,
																			event->ph[i].etas2(),
																			event->ph[i].etap(),
																			event->ph[i].cl_eta(),
    	                                                                	40,
    	                                                                	isMC,
																			event->ph[i].topoEtcone40(),
    	                                                                	false, 
    	                                                                	CaloIsoCorrection::PHOTON,
    	                                                                	CaloIsoCorrection::REL17_2);
			ph_Etcone40_corr->push_back(caloIsoCorr);			
		}

	}

	// cluster pt for 'electron' fsr
	vector<float>* el_bfEP_cl_pt = new vector<float>();	
	// Sanity check
	if( el_cur->n() != (Int_t) corr->bfEP_cl_pt.size() && doCorr) cout<<"Error: beforeEP_cl_pt has a different size than electron container"<<endl;
	for(Int_t i = 0; i < el_cur->n(); i++)
	{
		if(doCorr) el_bfEP_cl_pt->push_back(corr->bfEP_cl_pt[i]);
		else el_bfEP_cl_pt->push_back((*el_cur)[i].cl_pt());
	}

	// looping over the leptons and getting the FSR photons
	vector<ChargedLepton *>::iterator itr_i;
	Int_t i = -1;
	for(itr_i = lep_cont.begin(); itr_i != lep_cont.end(); ++itr_i)
	{
		i++;
		ChargedLepton *lep_curr = *itr_i;

		//// Topo ET cone
		//vector<float>* ph_Etcone40_corrected;
		//if(dataYear == 2011) ph_Etcone40_corrected = event->ph.Etcone40_pt_corrected();
		//else if(dataYear == 2012) ph_Etcone40_corrected = event->ph.topoEtcone40_corrected();
		
		if(lep_curr->getFlavor() == flavor::Muon)
		{	
			D3PDReader::MuonD3PDObjectElement* mu_curr = lep_curr->GetMuon();
			FsrPhotons fsrphotons;
			FsrPhotons::FsrCandidate cand_curr;

			if(lep_curr->getType() == leptonType::MuonCalo) 
			{
				candVec.push_back(cand_curr);
				continue;
			}

			if (fsrphotons.getFsrPhotons(mu_curr->eta(), 
    		                             mu_curr->phi(), 
    		                             mu_curr->id_theta(), 
    		                             mu_curr->id_phi(), 
    		                             event->ph.eta(), 
    		                             event->ph.phi(), 
    		                             event->ph.Et(), 
    		                             event->ph.f1(),
										 ph_Etcone40_corr,
    		                             event->ph.author(), 
    		                             el_cur->cl_eta(), 
    		                             el_cur->cl_phi(), 
    		                     	     el_bfEP_cl_pt, 
    		                             el_cur->f1(), 
    		                             el_cur->tracktheta(), 
    		                             el_cur->trackphi(),
    		                             cand_curr,
										 event->ph.tight())){
				isDebugCall = true;
				if(isDebugCall) {cout<<"Muon "<< i <<" Muon eta: "<<mu_curr->eta()<<" Muon Phi"<<mu_curr->phi()<<endl;
    		      std::cout << "Found FSR index, dR, Et, f1, cont " << cand_curr.index << " " << cand_curr.deltaR
    		                 << " " << cand_curr.Et << " " <<  cand_curr.f1 << " " << cand_curr.container << std::endl;
				cout<<"-------------------"<<endl;}
				isDebugCall = false;
    		  }
			candVec.push_back(cand_curr);
		}
		else if(lep_curr->getFlavor() == flavor::Electron)
		{	
			D3PDReader::ElectronD3PDObjectElement* el_curr = lep_curr->GetElectron();
			FsrPhotons fsrphotons;
			FsrPhotons::FsrCandidate cand_curr;
			if (fsrphotons.getFsrPhotons(el_curr->tracketa(),
										 el_curr->trackphi(),
    		                             event->ph.eta(), 
    		                             event->ph.phi(), 
    		                             event->ph.Et(), 
    		                             event->ph.f1(),
										 ph_Etcone40_corr,										   
    		                             event->ph.author(), 
    		                             event->ph.tight(), 
      		                             cand_curr)){
				if(isDebugCall) {cout<<"electron "<< i <<endl;
    		      std::cout << "Found FSR index, dR, Et, f1, cont " << cand_curr.index << " " << cand_curr.deltaR
    		                 << " " << cand_curr.Et << " " <<  cand_curr.f1 << " " << cand_curr.container << std::endl;
				cout<<"-------------------"<<endl;}
    		  }
			candVec.push_back(cand_curr);	
		}
	}

	// clearing up memory
	el_bfEP_cl_pt->clear();
	ph_Etcone40_corr->clear();
	delete ph_Etcone40_corr;
	delete el_bfEP_cl_pt;

	// Sanity check
	if(candVec.size() != 4) {cout<< "Correct FSR error: Don't have 4 FSR candidate"<<endl;}

	// Have found all the fsr Candidates
	// Do the analysis first for mumu Z1
	Bool_t passFSRCollinear = false;
	Bool_t passFSRFar = false;
	
	// Collinear only happens for 4mu or 2mu2el
	if(Type == analysisType::Mu4 || Type == analysisType::Mu2El2)
	{
		Double_t higestEt = -99999;
		FsrPhotons::FsrCandidate candidate;
		Bool_t isFSR = false;
		// Look at just the first two FSR and find the highest cluster Et fsr
		for(Int_t i = 0; i < 2; i++)
		{
			FsrPhotons::FsrCandidate tempCandidate = candVec[i];

			// Find the higest Et collinear photon on the cluster
			if(tempCandidate.index != -1 && tempCandidate.Et > higestEt && tempCandidate.deltaR <= 0.15)
			{
				higestEt = tempCandidate.Et;
				candidate = tempCandidate;
				isFSR = true;
			}
		}
		isDebugCall = true;
		if(isDebugCall) cout<<"FSR colinear candidate Passed the deltaR cut"<<endl;
		// if FSR found, then do the checks associated with Z1 for the collinear FSR
		if(isFSR)
		{
			if((higgs->getZ1()->get4Momentum()->M() < 89*1000) && (higgs->getZ1()->get4Momentum()->M() > 66*1000))
			{	
				if(isDebugCall) cout<<"FSR collinear candidate passed the Z1 mass cut"<<endl;
    		    TLorentzVector momFSR;
				momFSR.SetPtEtaPhiM(candidate.Et,candidate.eta,candidate.phi, 0.0);
	

				Double_t combMass = (higgs->getZ1()->getLepPlus()->get4MomentumNoP() 
						+ higgs->getZ1()->getLepNeg()->get4MomentumNoP() 
						+ momFSR).M();
				if(combMass < 100*1000)
				{
					if(isDebugCall) cout<<"FSR collinear candidate passed the combined mass cut"<<endl;	

					PATCore::ParticleType::Type particleType = PATCore::ParticleType::Electron; 
					if(candidate.container == "electron") particleType = PATCore::ParticleType::Electron;
					else if(candidate.container == "photon")
					{
						if(event->ph[candidate.index].isConv()) particleType = PATCore::ParticleType::ConvertedPhoton;
						else particleType = PATCore::ParticleType::UnconvertedPhoton;
					}
					higgs->FillFSRCorr(momFSR, runNumber_sf, true, false, particleType);
					passFSRCollinear = true;
					higgs->fsrType = fsrType::collFSRZ1mumu;				
				}
			}
		}
		isDebugCall = false;
	}
	// If no collinear FSR has been found
	if(!passFSRCollinear)
	{
		
		Double_t higestEt = -99999;
		FsrPhotons::FsrCandidate candidate;
		Bool_t isFSR = false;
		// Loop over all the candidates and find the highest Et cluster far FSR
		for(Int_t i = 0; i < (Int_t) candVec.size(); i++)
		{
			FsrPhotons::FsrCandidate tempCandidate = candVec[i];

			// Find the higest Et far photon on the cluster
			if(tempCandidate.index != -1 && tempCandidate.Et > higestEt && tempCandidate.deltaR > 0.15)
			{
				higestEt = tempCandidate.Et;
				candidate = tempCandidate;
				isFSR = true;
			}
		}
		if(isDebugCall) cout<<"FSR far candidate exists"<<endl;					
		
		// If there is such a FSR candidate
		if(isFSR)
		{
				//cout<<"There is a candidate Far FSR"<<endl;
				TLorentzVector momFSR;
				momFSR.SetPtEtaPhiM(candidate.Et,candidate.eta,candidate.phi, 0.0);

				// Delta R cut: Ensure that deltaR > 0.15 for all the leptons in the quadruplet
				// Sortof the opposite logic here
				Bool_t deltaRCut = false;
				for(itr_i = lep_cont.begin(); itr_i != lep_cont.end(); ++itr_i)
				{
					ChargedLepton *lep_curr = *itr_i;
					//cout<<"delta R: "<<momFSR.DeltaR(lep_curr->get4MomentumNoP())<<endl;
					if(momFSR.DeltaR(lep_curr->get4MomentumNoP()) < 0.15) deltaRCut = true;
				}
				// If deltaR cut has been pases
				if(!deltaRCut)
				{
					if(isDebugCall) cout<<"FSR far candidate passed Delta R cut"<<endl;										
					//cout<<"Passed Delta R cut"<<endl;
					Double_t deltaZ1Mass = 999999;
					Double_t deltaZ2Mass = 999999;
					// if the dilepton mass is less than 81 GeV, find how much the FSR move the mass close to 
					// Z pole mass... 100 GeV limit there as well
					if(higgs->getZ1()->get4Momentum()->M() < 81*1000) 
					{
						TLorentzVector combined = higgs->getZ1()->get4MomentumNoP() + momFSR;
						if(combined.M() < 100 * 1000) deltaZ1Mass = fabs(combined.M() - pdgZMass);
					}
					if(higgs->getZ2()->get4Momentum()->M() < 81*1000) 
					{
						TLorentzVector combined = higgs->getZ2()->get4MomentumNoP() + momFSR;
						if(combined.M() < 100 * 1000) deltaZ2Mass = fabs(combined.M() - pdgZMass);
						if(higgs->get4Momentum()->M() < 190*1000) deltaZ2Mass = 999999;
					}
					if(isDebugCall) cout<<"deltaZ1Mass: "<<deltaZ1Mass<<" deltaZ2Mass: "<<deltaZ2Mass<<endl;
					// Figuring out which dilepton to associate the far FSR photon with
					if(deltaZ1Mass < deltaZ2Mass && deltaZ1Mass != 999999) 
					{
						PATCore::ParticleType::Type particleType = PATCore::ParticleType::Electron; 
						if(candidate.container == "electron") particleType = PATCore::ParticleType::Electron;
						else if(candidate.container == "photon")
						{
							if(event->ph[candidate.index].isConv()) particleType = PATCore::ParticleType::ConvertedPhoton;
							else particleType = PATCore::ParticleType::UnconvertedPhoton;
						}
						higgs->FillFSRCorr(momFSR, runNumber_sf, true, false, particleType);
						passFSRFar = true;
						if(isDebugCall) cout<<"Far FSR Attached to Z1"<<endl;
						higgs->fsrType = fsrType::farFSRZ1;
					}
					else if (deltaZ2Mass < deltaZ1Mass && deltaZ2Mass != 999999) 
					{
						PATCore::ParticleType::Type particleType = PATCore::ParticleType::Electron; 
						if(candidate.container == "electron") particleType = PATCore::ParticleType::Electron;
						else if(candidate.container == "photon")
						{
							if(event->ph[candidate.index].isConv()) particleType = PATCore::ParticleType::ConvertedPhoton;
							else particleType = PATCore::ParticleType::UnconvertedPhoton;
						}
						higgs->FillFSRCorr(momFSR, runNumber_sf, false, true, particleType);
						passFSRFar = true;						
						if(isDebugCall) cout<<"Far FSR attached to Z2"<<endl;
						higgs->fsrType = fsrType::farFSRZ2;
					}
				}
		}
	}

	if(passFSRCollinear && isDebugCall)cout<<"-------------------"<<endl<<"Colinear FSR event"<<endl<<"-------------------"<<endl;
	if(passFSRFar && isDebugCall)cout<<"-------------------"<<endl<<"Far FSR event"<<endl<<"-------------------"<<endl;
	// Calculating the FSR mass...
	Double_t mass = InvariantMass::invMassCalc(higgs->getLeptonFSRLorentz());
	higgs->SetMassFSR(mass);
	Double_t massErr = MassError::sigmaMassCalc_d0z0PhiThetaP(higgs->getLeptonFSRLorentz(),higgs->getLeptonFSRCovMatrix());
	higgs->SetMassErrFSR(massErr);
	
	// Z1 Mass
	Double_t z1Mass = higgs->getZ1()->get4Momentum()->M();
	if(higgs->getZ1()->getHasFSR())
	{
		z1Mass = (higgs->getZ1()->get4MomentumNoP() + higgs->getFSR4MomentumNoP()).M();
	}
	higgs->SetZ1FSRMass(z1Mass);

	// Z2 mass
	Double_t z2Mass = higgs->getZ2()->get4Momentum()->M();
	if(higgs->getZ2()->getHasFSR())
	{
		z2Mass = (higgs->getZ2()->get4MomentumNoP() + higgs->getFSR4MomentumNoP()).M();
		//cout<<event->eventinfo.EventNumber()<<" Correct Mass: "<<z2Mass/1000<<" old mass: "<<higgs->getZ2()->get4MomentumNoP().M()/1000<<endl;
	}
	higgs->SetZ2FSRMass(z2Mass);

	// Creating a sum of TLorentzVector for the overall sum
	vector<TLorentzVector> leptonLorentzFsr = higgs->getLeptonFSRLorentz();
	for(Int_t i = 0; i <  (Int_t) leptonLorentzFsr.size(); i++)
	{
		higgs->sum_fsr = higgs->sum_fsr + leptonLorentzFsr[i];
	}
	leptonLorentzFsr.clear();
	
	isDebugCall = false;
}	
// Zmass Constraint
void HiggsAnalysis::CorrectZMassConstraint(QuadLepton * higgs, Int_t muType)
{
	Bool_t hasWidth = true;
	ZMassConstraint::ConstraintFit *massFit=new ZMassConstraint::ConstraintFit(pdgZMass,hasWidth,2.4952e3);
	// For Stroage
	vector <TLorentzVector> lep4Momentum;
	vector <TMatrixD> lepErr;
	// for Z1
	vector <TLorentzVector> Z14Momentum;
	// for Z1
	vector <TLorentzVector> Z24Momentum;

	// Zmass Constraint for Z1
    ZMassConstraint::ConstraintFitInput inputZMassConstraint;
	inputZMassConstraint.addConstituent_FourVector_d0z0PhiThetaP(higgs->getZ1()->getLepPlus()->get4MomentumNoP(muType),
                                                           higgs->getZ1()->getLepPlus()->getHepCovMatrix(muType));
          
    inputZMassConstraint.addConstituent_FourVector_d0z0PhiThetaP(higgs->getZ1()->getLepNeg()->get4MomentumNoP(muType),
                                                           higgs->getZ1()->getLepNeg()->getHepCovMatrix(muType));

	// Add the FSR Photon it is there
	if(higgs->getZ1()->getHasFSR())
	{
		inputZMassConstraint.addConstituent_FourVector_d0z0PhiThetaP(higgs->getZ1()->getFSR4Momentum(),
                                                           higgs->getZ1()->getHepFSRError());
	}

	// Perform the Zmass contraint on Z1
	massFit->massFitInterface(inputZMassConstraint);
    ZMassConstraint::ConstraintFitOutput massFitResult = massFit->massFitRun(-1.);
	
	// Storing the Results
	lep4Momentum.push_back(massFitResult.getConstituentFourVector(0));
	lep4Momentum.push_back(massFitResult.getConstituentFourVector(1));
	lepErr.push_back(massFitResult.getConstituentCovarianced0z0PhiThetaP_TMatrix(0));
	lepErr.push_back(massFitResult.getConstituentCovarianced0z0PhiThetaP_TMatrix(1));
	// Getting the FSR Photon
	if(higgs->getZ1()->getHasFSR() && massFitResult.getNConstituents()>2)
	{
		lep4Momentum.push_back(massFitResult.getConstituentFourVector(2));
		lepErr.push_back(massFitResult.getConstituentCovarianced0z0PhiThetaP_TMatrix(2));
	} // Sanity check
	else if (higgs->getZ1()->getHasFSR() && !(massFitResult.getNConstituents()>2))
	{
		cout<<"ERROR: CorrectZMassConstraint: Z1 FSR Event but no Photon correction from ZMass Constraints"<<endl;
		cout<<"Has FSR: "<< higgs->getZ1()->getHasFSR()<< " n Zmass vector: "<<massFitResult.getNConstituents()<<endl<<endl;
	}
	Z14Momentum = lep4Momentum;
		
	// Do Z mass constraint on Z2 if the mass is above 190 Gev
	if(higgs->get4Momentum()->M() > 190*1000)
	{
    	ZMassConstraint::ConstraintFitInput inputZ2MassConstraint;
		inputZ2MassConstraint.addConstituent_FourVector_d0z0PhiThetaP(higgs->getZ2()->getLepPlus()->get4MomentumNoP(muType),
                                                           higgs->getZ2()->getLepPlus()->getHepCovMatrix(muType));
          
    	inputZ2MassConstraint.addConstituent_FourVector_d0z0PhiThetaP(higgs->getZ2()->getLepNeg()->get4MomentumNoP(muType),
                                                           higgs->getZ2()->getLepNeg()->getHepCovMatrix(muType));

		// Add the FSR Photon it is there
		if(higgs->getZ2()->getHasFSR())
		{
			inputZ2MassConstraint.addConstituent_FourVector_d0z0PhiThetaP(higgs->getZ2()->getFSR4Momentum(),
    	                                                       higgs->getZ2()->getHepFSRError());
		}

		// Perform the Zmass contraint on Z2
		massFit->massFitInterface(inputZ2MassConstraint);
    	ZMassConstraint::ConstraintFitOutput massFitZ2Result = massFit->massFitRun(-1.);
		
		// Storing the Results
		lep4Momentum.push_back(massFitZ2Result.getConstituentFourVector(0));
		lep4Momentum.push_back(massFitZ2Result.getConstituentFourVector(1));
		lepErr.push_back(massFitZ2Result.getConstituentCovarianced0z0PhiThetaP_TMatrix(0));
		lepErr.push_back(massFitZ2Result.getConstituentCovarianced0z0PhiThetaP_TMatrix(1));
		Z24Momentum.push_back(massFitZ2Result.getConstituentFourVector(0));
		Z24Momentum.push_back(massFitZ2Result.getConstituentFourVector(1));

		// Getting the FSR Photon
		if(higgs->getZ2()->getHasFSR() && massFitZ2Result.getNConstituents()>2)
		{
			lep4Momentum.push_back(massFitZ2Result.getConstituentFourVector(2));
			lepErr.push_back(massFitZ2Result.getConstituentCovarianced0z0PhiThetaP_TMatrix(2));
			Z24Momentum.push_back(massFitZ2Result.getConstituentFourVector(2));	
		} // Sanity check
		else if (higgs->getZ2()->getHasFSR() && !(massFitZ2Result.getNConstituents()>2))
		{
			cout<<"ERROR: CorrectZMassConstraint: Z2 FSR Event but no Photon correction from ZMass Constraints"<<endl;
			cout<<"Has FSR: "<< higgs->getZ2()->getHasFSR()<< " n Zmass vector: "<<massFitResult.getNConstituents()<<endl<<endl;
		}
	}
	else
	{
		// If no Z mass constraints for Z2, add just the normal 4Vector and Cov Matrix 
		lep4Momentum.push_back(higgs->getZ2()->getLepPlus()->get4MomentumNoP(muType));
		lep4Momentum.push_back(higgs->getZ2()->getLepNeg()->get4MomentumNoP(muType));
		lepErr.push_back(higgs->getZ2()->getLepPlus()->getCovMatrix(muType));
		lepErr.push_back(higgs->getZ2()->getLepNeg()->getCovMatrix(muType));
		Z24Momentum.push_back(higgs->getZ2()->getLepPlus()->get4MomentumNoP(muType));
		Z24Momentum.push_back(higgs->getZ2()->getLepNeg()->get4MomentumNoP(muType));

		// Add the fsr if it is there as well
		if(higgs->getZ2()->getHasFSR())
		{
			lep4Momentum.push_back(higgs->getZ2()->getFSR4Momentum());
			lepErr.push_back(higgs->getZ2()->getFSRError());
			Z24Momentum.push_back(higgs->getZ2()->getFSR4Momentum());
		}
	}

	// Stroing the event in the Higgs event
	if(muType == muonType::CB) higgs->SetZmassLorentz(lep4Momentum);
	if(muType == muonType::CB) higgs->SetZmassCovMat(lepErr);

	// Getting the mass and storing it
	//Double_t mass = InvariantMass::invMassCalc(higgs->getLeptonZMassLorentz());
	Double_t mass = InvariantMass::invMassCalc(lep4Momentum);
	higgs->SetMassZMass(mass, muType);
	//Double_t massErr = MassError::sigmaMassCalc_d0z0PhiThetaP(higgs->getLeptonZMassLorentz(),higgs->getLeptonZMassCovMatrix());
	Double_t massErr = MassError::sigmaMassCalc_d0z0PhiThetaP(lep4Momentum,lepErr);
	higgs->SetMassErrZMass(massErr, muType);
	Double_t z1Mass = InvariantMass::invMassCalc(Z14Momentum);
	higgs->SetZ1MassZMass(z1Mass, muType);
	Double_t z2Mass = InvariantMass::invMassCalc(Z24Momentum);
	higgs->SetZ2MassZMass(z2Mass, muType);

	// Creating a sum of TLorentzVector for the overall sum
	vector<TLorentzVector> leptonLorentzZmass = higgs->getLeptonZMassLorentz();
	for(Int_t i = 0; i <  (Int_t) leptonLorentzZmass.size(); i++)
	{
		if(muType == muonType::CB) higgs->sum_constrained = higgs->sum_constrained + leptonLorentzZmass[i];
	}

	lep4Momentum.clear();
	lepErr.clear();
	Z14Momentum.clear();
	Z24Momentum.clear();
	leptonLorentzZmass.clear();

	delete massFit;
}

// Gets the event weight
// Zvertex, pileup and overlap with ttbar, bb, and ZZ
Double_t HiggsAnalysis::getEventWeight()
{
	Double_t weight = 1;
	if(doZVertexWeight) weight = weight * getZVertexWeight();
	if(doPileupWeight) weight = weight * getPileupWeight();
	if(doggFWeight) weight = weight * getggFWeight();
	if(doJHUWeight) weight = weight * getJHUWeight();

	//if(doZZOverlap) weight = weight * getZZOverlapVeto();
	//if(doTTBarOverlap) weight = weight * getTTBarVeto();
	//if(doZBBOverlap) weight = weight * getZBBOverlapVeto();
	// Adding the MC weight
	weight = weight * event->eventinfo.mc_event_weight();
	return weight;
}
// ZVertex weight
Double_t HiggsAnalysis::getZVertexWeight()
{
	Double_t zWeight = 1;
	// Protecting against Data
	if(!isMC) return 1;
	// Protection for when combining things together
	if(!doZVertexWeight) return zWeight;

	if(event->mc.n() > 1)
	{
		Double_t vxp_z = event->mc[2].vx_z();
		zWeight = zVertexTool->GetWeight(vxp_z);
	}
	return zWeight;
}
// ggF reweighting
Double_t HiggsAnalysis::getggFWeight()
{
	Double_t ggFWeight = 1;
	// Protection for when combining things together
	if(!doggFWeight || !isMC || ggFReweight == 0) return ggFWeight;
	// To get the true Higgs Pt
	//Int_t higgsPdgId = 25;
	//Int_t dataset = event->eventinfo.mc_channel_number();
	//// Dependant on data
	//if (dataset == 167120 || dataset == 167121 || dataset == 167122 || dataset == 167123) {
	//   higgsPdgId = 39;
	//
	Double_t trueHiggsPt = 0;
	
   	for (Int_t i = 0; i < event->mc.n(); i++) {
		Int_t pdgID = TMath::Abs(event->mc[i].pdgId());
		Int_t status = event->mc[i].status();
		
		//if ( ((pdgID == 25 || pdgID == 39) ) &&(status==2 || status==10902 || status==62)) {
		if ( ((pdgID == 25 ) ) &&(status==2 || status==10902 || status==62)) {
			trueHiggsPt = event->mc[i].pt()/1000;
		}
	}
	//if(isDebugCall)
	//{
	//	cout<<"-----------------"<<endl;
	//	cout<<"ggF reweight"<<endl;
	//	cout<<"Higgs Pt: "<<trueHiggsPt<<endl;
	//	cout<<"-----------------"<<endl;

	//}
	pair<double, double> result;

	if(sampleProdType == sampleType::ggF || sampleProdType == sampleType::ggF_ZpZp) 
	{
		result = ggFReweight->getWeightAndStatError(trueHiggsPt);
		ggFWeight=result.first;
	}

	return ggFWeight;
	

}
// JHU reweighting
Double_t HiggsAnalysis::getJHUWeight()
{
	Double_t JHUWeight = 1;
	// Protection for when combining things together
	if(!doJHUWeight) return JHUWeight;
	// To get the true Higgs Pt
	Double_t trueHiggsPt = 0;
	
   		for (Int_t i = 0; i < event->mc.n(); i++) {
		Int_t pdgID = TMath::Abs(event->mc[i].pdgId());
		Int_t status = event->mc[i].status();
		
		if ( ((pdgID == 25) || (pdgID == 39)) &&(status==1 || status==10902 || status==62)) {
			trueHiggsPt = event->mc[i].pt();
		}
	}

	JHUWeight = ptJHUReweighting->GetJHUWeight(event->eventinfo.mc_channel_number(), trueHiggsPt/1000, 0);
	
	return JHUWeight;
	

}

// PileUp Weight
Double_t HiggsAnalysis::getPileupWeight()
{
	// Average crossing per interaction
	// Protecting against Data
	if(!isMC) return 1;

	Double_t mu = -9999;
	mu = event->eventinfo.averageIntPerXing();
	if(isMC)
	{
		if(dataYear == 2012 && isMC)  mu=(event->eventinfo.lbn()==1&&int(mu+0.5)==1)?0.:mu;
	}
	// mu=round(mu);
	// Just a  place holder
	Double_t pileupWeight = 1;
	// Getting the actual weight
	pileupWeight = pileupTool->GetCombinedWeight(event->eventinfo.RunNumber(),event->eventinfo.mc_channel_number(),mu);

	return pileupWeight;;
}
// Veto for TTbar Overlap
Double_t HiggsAnalysis::getTTBarVeto()
{
	// Protecting against Data
	if(!isMC) return 1;

	Double_t veto = 1;

	Int_t mcChannelNum1 = event->eventinfo.mc_channel_number();
	// mc_channel_number2 : core sample                        
    // 146369 to set up the 2012 ttbar filter
    // 109346 to set up the 2011 ttbar filter
	Int_t mcChannelNum2 = 111111;
	if(isMC)
	{
		if(dataYear == 2011) mcChannelNum2 = 109346;
		else if(dataYear == 2012) mcChannelNum2 = 146369;
	}
	// Function from HiggsZZ4Utils
	Bool_t ttbarOverlap = killEvent(mcChannelNum1, 
									mcChannelNum2,
									event->mc.n(),
									event->mc.pt(),
									event->mc.eta(),
									event->mc.phi(),
									event->mc.m(),
									event->mc.status(),
									event->mc.pdgId());
	if(ttbarOverlap) veto = 0;

	return veto;
}
// Veto for TTbar Overlap
Double_t HiggsAnalysis::getZZOverlapVeto()
{
	// Protecting against Data
	if(!isMC) return 1;

	Double_t veto = 1;

	Int_t mcChannelNum1 = event->eventinfo.mc_channel_number();
	Int_t mcChannelNum2 = 111111;
	// Function from HiggsZZ4Utils
	Bool_t zzOverlap = killEvent(mcChannelNum1, 
								mcChannelNum2,
								event->mc.n(),
								event->mc.pt(),
								event->mc.eta(),
								event->mc.phi(),
								event->mc.m(),
								event->mc.status(),
								event->mc.pdgId());
	if(zzOverlap) veto = 0;

	return veto;
}
// ZBB Overlap
Double_t HiggsAnalysis::getZBBOverlapVeto()
{
	// Protecting against Data
	if(!isMC) return 1;

	Double_t veto = 1;

	if (event->top.hfor_type() == 4)
	{
        veto = 0.;
    }
	return veto;
}

// Weights specific to the higgs Candidate
Double_t HiggsAnalysis::getHiggsWeight(QuadLepton * higgs)
{
	Double_t weight = 1;

	weight = weight * getTriggerWeight(higgs);
	weight = weight * getLepEffWeight(higgs);

	return weight;
}

// Weights that come from single trigger efficiency
Double_t HiggsAnalysis::getTriggerWeight(QuadLepton * higgs)
{
	Double_t trigWeight = 1;

	Int_t analysisType = higgs->getQuadType();
	Bool_t singleTriggered = false;
	// Checking if the single trigger was triggered
	if(analysisType == quadType::Mu4) {singleTriggered = SingleMuonTrigger(runNumber_sf);}
	else if (analysisType == quadType::El4) {singleTriggered = SingleElectronTrigger(runNumber_sf);}
	else if (analysisType == quadType::El2Mu2 || analysisType == quadType::Mu2El2)
	{
		singleTriggered = SingleMuonTrigger(runNumber_sf) || SingleElectronTrigger(runNumber_sf);
	}
	else {cout<<"Error: HiggsAnalysis::getTriggerWeight: Quad Type not Recognized";}

	if(singleTriggered)
	{
		// Variable for trigger effciency
        vector<TLorentzVector>    elTrigSF;
        vector<TLorentzVector>    muTrigSF;
        vector<muon_quality>      muQuality;
        vector<electron_quality>  elQuality;
		
		// Setting Electron quality
		electron_quality p;
		if(electronCollection == electronCollection::LoosePlusPlus) p = loosepp;
		else if(electronCollection == electronCollection::MultiLepton) p = ML;
		else if (electronCollection == electronCollection::Likelihood) {p = LooseLLH; }//cout<<"Trigger SF selection needs to be fixed"<<endl;}
		else {cout<<"Error: HiggsAnalysis::getTriggerWeight: Electron collection not recognized"<<endl;}
		
		// Muon quality
		muon_quality q;
		if(muonCollection == muonCollection::Loose) q = loose;
		else {cout<<"Error: HiggsAnalysis::getTriggerWeight: Muon collection not recognized"<<endl;}
		
		// Filling the vectors
		vector<ChargedLepton *> lep = higgs->getLepton();
		for(vector<ChargedLepton *>::iterator itr_lep = lep.begin();
			itr_lep != lep.end(); ++itr_lep)
		{
			ChargedLepton* curr_i = *itr_lep;
			if(curr_i->getFlavor() == flavor::Electron)
			{
				elTrigSF.push_back(curr_i->get4MomentumTrigSF());
				elQuality.push_back(p);
			}
			else if(curr_i->getFlavor() == flavor::Muon)
			{
				muTrigSF.push_back(curr_i->get4MomentumTrigSF());
				if(curr_i->getType() != leptonType::MuonCalo) muQuality.push_back(q);
				else if (curr_i->getType() == leptonType::MuonCalo) muQuality.push_back(CaloMuon);
			}
		}
		
		// Getting the trigger Efficiency
		pair<double,double> SF = leptonSF->GetTriggerSF(runNumber_sf,
                                                        false, 
                                                        muTrigSF, 
                                                        muQuality,    //q, 
                                                        elTrigSF, 
                                                        elQuality, //p, 
                                                        0);

		if(SF.first != 0.0) {trigWeight = SF.first;}
		else {cout<<"HiggsAnalysis::getTriggerWeight: Problem with Trigger SF"<<endl;}
      	elTrigSF.clear();
        muTrigSF.clear();
        muQuality.clear();
        elQuality.clear();
		lep.clear();
	}

	return trigWeight;
}
// ID and Reco Weight
Double_t HiggsAnalysis::getLepEffWeight(QuadLepton * higgs)
{
	Double_t lepWeight = 1;
	vector<ChargedLepton *> lep = higgs->getLepton();
	for(vector<ChargedLepton *>::iterator itr_lep = lep.begin();
		itr_lep != lep.end(); ++itr_lep)
	{
		ChargedLepton* curr_i = *itr_lep;
		lepWeight = lepWeight * curr_i->getLepEff();
	}
	
	lep.clear();
	return lepWeight;
}
// Filling the weights
void HiggsAnalysis::fillHiggsWeight(QuadLepton * higgs)
{
	if(isMC)
	{
		higgs->zVertWeight = getZVertexWeight();
		higgs->pileupWeight = getPileupWeight();
		higgs->ggFWeight = getggFWeight();		
		higgs->JHUWeight = getJHUWeight();		
		higgs->ttBarVeto = getTTBarVeto();
		higgs->zzOverlapVeto = getZZOverlapVeto();
		higgs->zBBVeto = getZBBOverlapVeto();
		higgs->trigEff = getTriggerWeight(higgs);
		higgs->lepEff = getLepEffWeight(higgs);
		higgs->mcEventWeight = event->eventinfo.mc_event_weight();
	}
	else
	{
		higgs->zVertWeight = 1;
		higgs->pileupWeight = 1;
		higgs->ggFWeight = 1;				
		higgs->JHUWeight = 1;				
		higgs->ttBarVeto = 1;
		higgs->zzOverlapVeto = 1;
		higgs->zBBVeto = 1;
		higgs->trigEff = 1;
		higgs->lepEff = 1;
		higgs->mcEventWeight = 1;
	}

	if(printWeight || isDebugCall)
	{
		Int_t printL = 25;
		Int_t pres = 10;

		Double_t zVertexWeight = 1;
		Double_t ggFWeight = 1;
		Double_t JHUWeight = higgs->JHUWeight;
		
		if(doZVertexWeight)zVertexWeight = getZVertexWeight();
		if(doggFWeight)ggFWeight = getggFWeight();
		
		cout<<setw(printL)<<"pu weight: "<<setprecision(pres)<<getPileupWeight()<<endl;
		cout<<setw(printL)<<"vxz weight: "<<setprecision(pres)<<zVertexWeight<<endl;
		cout<<setw(printL)<<"trigSF weight: "<<setprecision(pres)<<getTriggerWeight(higgs)<<endl;
		cout<<setw(printL)<<"ggF weight: "<<setprecision(pres)<<ggFWeight<<endl;
		cout<<setw(printL)<<"JHU weight: "<<setprecision(pres)<<JHUWeight<<endl;
		cout<<setw(printL)<<"mc Event weight: "<<setprecision(pres)<<higgs->mcEventWeight<<endl;
		vector<ChargedLepton *> lep = higgs->getLepton();
		Double_t lepWeight = 1;
		for(vector<ChargedLepton *>::iterator itr_lep = lep.begin();
			itr_lep != lep.end(); ++itr_lep)
		{
			ChargedLepton* curr_i = *itr_lep;
			lepWeight = lepWeight * curr_i->getLepEff();
			TString lepType = "";
			if(curr_i->getFlavor() == flavor::Muon) lepType = "Muon ";
			else if (curr_i->getFlavor() == flavor::Electron) lepType = "Electron ";
			lepType = lepType + "Lep eff weight: ";
			cout<<setw(printL)<<lepType<<setprecision(pres)<<curr_i->getLepEff()<<endl;
		}
		cout<<setw(printL)<<"lep Eff: "<<setprecision(pres)<<lepWeight<<endl;
		cout<<setw(printL)<<"TotalWeight: "<<setprecision(pres)<<getPileupWeight()*zVertexWeight*getTriggerWeight(higgs)*lepWeight*ggFWeight*JHUWeight<<endl;
		cout<<setw(printL)<<"ttbar veto: "<<setprecision(pres)<<higgs->ttBarVeto<<endl;
		cout<<setw(printL)<<"zzOverlapVeto: "<<setprecision(pres)<<higgs->zzOverlapVeto<<endl;
		cout<<setw(printL)<<"zBB veto: "<<setprecision(pres)<<higgs->zBBVeto<<endl;
	
	}
}

// CrossSection Stuff
void HiggsAnalysis::fillCrossSection(QuadLepton * higgs)
{
	// Setting up the vars
	Double_t crossSection = -1;
	Double_t BR_correction_4l = -1;
	Double_t BR_correction_2l2l = -1;
	Double_t BR_correction = -1;
	Double_t lumi = -1;

	if(isMC)
	{
		// Getting the run number
		Int_t RunNumber  = event->eventinfo.mc_channel_number();
		if(mcRunNumber != RunNumber)
		{
			cout<<"Error: getCrossSection: D3PD runNumber("<<RunNumber<<") doesn't match the file runNumber ("<<mcRunNumber<<")"<<endl;
		}
		Double_t massHiggs = getMCHiggsMass();
		CrossSections::LHCEnergy energyLHC;
		if(dataYear == 2011) energyLHC = CrossSections::SevenTeV;
		else if(dataYear == 2012) energyLHC = CrossSections::EightTeV;
	
		// get cross section 
		if( (RunNumber>=116761 && RunNumber<=116799) || // 2011
            (RunNumber>=116611 && RunNumber<=116612) || // 2011
            (RunNumber>=116994 && RunNumber<=116997) || // 2011
            (RunNumber>=146270 && RunNumber<=146275) || // 2011
            (RunNumber>=169494 && RunNumber<=169497) || // 2011
            (RunNumber>=167600 && RunNumber<=167603) || // 2011 JHU
            (RunNumber>=167606 && RunNumber<=167606) || // 2011 JHU
           
            (RunNumber>=160351 && RunNumber<=160354) || // 2012 Additional
            (RunNumber>=167876 && RunNumber<=167878) || // 2012 Additional
            (RunNumber>=160180 && RunNumber<=160190) || // 2012 Additional
            (RunNumber>=160155 && RunNumber<=160155) || // 2012 Additional
            (RunNumber>=167974 && RunNumber<=167974) || // 2012 Additional
 
            (RunNumber>=167120 && RunNumber<=167123) || // 2012 JHU
            (RunNumber>=167126 && RunNumber<=167126) || // 2012 JHU
            (RunNumber>=160152 && RunNumber<=160179) || // 2012
            (RunNumber>=167890 && RunNumber<=167894) || // 2012
            (RunNumber>=167220 && RunNumber<=167227) || // 2012            
            (RunNumber>=167960 && RunNumber<=167987) || // 2012
            (RunNumber>=181040 && RunNumber<=181055) || // 2012
			(RunNumber>=181330 && RunNumber<=181336) || // 2012 new samples
			(RunNumber>=167895 && RunNumber<=167896) || // 2012 new samples// added Oct 6
			(RunNumber>=181306 && RunNumber<=181307) || // 2012 new samples // added Oct 6
			(RunNumber>=181870 && RunNumber<=181874) || // 2012 new samples added Jan 3 14
			(RunNumber == 189140) 					 ){ // 2012 inclusive sample added Jan 31
			if(sampleProdType != sampleType::ggF) cout<<"Error: getCrossSection: prodType misMatch ggF"<<endl;	
			crossSection = Higgs_xs.higgsprodxsecGGF(massHiggs,energyLHC);
		}
		else if( (RunNumber>=125063 && RunNumber<=125101) ||  // 2011
                 (RunNumber>=116621 && RunNumber<=116622) ||  // 2011
                 (RunNumber>=169070 && RunNumber<=169073) ||  // 2011
                 (RunNumber>=146280 && RunNumber<=146285) ||  // 2011
                 (RunNumber>=169870 && RunNumber<=169873) ||  // 2011 
                 
				 (RunNumber>=160230 && RunNumber<=160230) || // 2012 Additional
                 (RunNumber>=160205 && RunNumber<=160205) || // 2012 Additional

				 (RunNumber>=160202 && RunNumber<=160229) ||  // 2012
                 (RunNumber>=160231 && RunNumber<=160240) ||  // 2012
                 (RunNumber>=167990 && RunNumber<=167996) ||  // 2012
                 (RunNumber>=167230 && RunNumber<=167237) ||  // 2012
                 (RunNumber>=169160 && RunNumber<=169187) ||  // 2012
                 (RunNumber>=181160 && RunNumber<=181175) ||  // 2012
                 (RunNumber>=181337 && RunNumber<=181343) ||  // 2012 new samples
				 (RunNumber>=181885 && RunNumber<=181889) ||  // 2012 samples added Jan 3 14
				 (RunNumber == 189141) 				 	  ){ // 2012 inclusive sample added Jan 31
			if(sampleProdType != sampleType::VBF) cout<<"Error: getCrossSection: prodType misMatch VBF"<<endl;			
		  	crossSection = Higgs_xs.higgsprodxsecVBF(massHiggs,energyLHC);
		} 
		else if( (RunNumber>=125265 && RunNumber<=125290) || // 2011
                 (RunNumber>=169874 && RunNumber<=169877) || // 2011
                 (RunNumber>=125800 && RunNumber<=125801) || // 2011
                 
                 (RunNumber>=160250 && RunNumber<=160280) || // 2012
                 (RunNumber>=167240 && RunNumber<=167247) || // 2012
				 (RunNumber == 189142) 				      ){ // 2012 inclusive sample added Jan 31
		  
			if(sampleProdType != sampleType::WH) cout<<"Error: getCrossSection: prodType misMatch WH"<<endl;	
			crossSection = Higgs_xs.higgsprodxsecWH(massHiggs,energyLHC);
		} 
		else if( (RunNumber>=125425 && RunNumber<=125450) || // 2011
                 (RunNumber>=125900 && RunNumber<=125901) || // 2011
                 (RunNumber>=169878 && RunNumber<=169881) || // 2011
                 
                 (RunNumber>=160300 && RunNumber<=160330) || // 2012
                 (RunNumber>=167250 && RunNumber<=167257) || // 2012
				 (RunNumber == 189143) 					  ){ // 2012 inclusive sample added Jan 31
		  
			if(sampleProdType != sampleType::ZH) cout<<"Error: getCrossSection: prodType misMatch ZH"<<endl;	
		  	crossSection = Higgs_xs.higgsprodxsecZH(massHiggs,energyLHC);
		}
		else if( (RunNumber>=167570 && RunNumber<=167575) || // 2011
 
                 (RunNumber>=167560 && RunNumber<=167569) || // 2012
				 (RunNumber == 189144) 					  ){ // 2012 inclusive sample added Jan 31
			
			if(sampleProdType != sampleType::ttH) cout<<"Error: getCrossSection: prodType misMatch ttH"<<endl;	
		  	crossSection = Higgs_xs.higgsprodxsecttH(massHiggs,energyLHC);
		}
		else if( (RunNumber>=169716 && RunNumber<=169717) ||  // 2011 JHU
                 (RunNumber>=167604 && RunNumber<=167605) ||  // 2011 JHU
                 (RunNumber>=167607 && RunNumber<=167607) ||  // 2012 JHU
  
                 (RunNumber>=169710 && RunNumber<=169711) ||  // 2012 JHU
                 (RunNumber>=167124 && RunNumber<=167125) ||  // 2012 JHU
                 (RunNumber>=167127 && RunNumber<=167127) ) {  // 2012 JHU		  
			if(sampleProdType != sampleType::qqF) cout<<"Error: getCrossSection: prodType misMatch qqF"<<endl;			
		  	crossSection = Higgs_xs.higgsprodxsecGGF(massHiggs,energyLHC);
		}		
		else if( (RunNumber>=181250 && RunNumber<=181280)) {  // 2012 z'z' s// added Oct 6
			if(sampleProdType != sampleType::ggF_ZpZp) cout<<"Error: getCrossSection: prodType misMatch  z'z' qqF"<<endl;	
			
			if(energyLHC == CrossSections::SevenTeV) crossSection = CrossSections::GetBkgCrossSection7TeV(RunNumber);
			else if (energyLHC == CrossSections::EightTeV) crossSection = CrossSections::GetBkgCrossSection8TeV(RunNumber);
		} 

		else{
		  //// HSG2 cross section
			if(sampleProdType != sampleType::Background) cout<<"Error: getCrossSection: prodType misMatch background"<<endl;	  
		  	//crossSection = bkgCrossSection(RunNumber,energyLHC);
			if(energyLHC == CrossSections::SevenTeV) crossSection = CrossSections::GetBkgCrossSection7TeV(RunNumber);
			else if (energyLHC == CrossSections::EightTeV) crossSection = CrossSections::GetBkgCrossSection8TeV(RunNumber);
		} 
		
		// No Tau sample
		if(noTauSample) 
		{
			//BR_correction_4l = higgs_bror->get4eBrRatio(massHiggs) / higgs_bror->getemutauBrRatio(massHiggs) * 4.;
	    	//BR_correction_2l2l = higgs_bror->get2e2muBrRatio(massHiggs) / higgs_bror->getemutauBrRatio(massHiggs) * 4. / 2.;
			BR_correction_4l = higgs_bror->get4eBrRatio(massHiggs)  * 4.;
	    	BR_correction_2l2l = higgs_bror->get2e2muBrRatio(massHiggs)  * 4. / 2.;
		}
		else 
		{
			//BR_correction_4l = higgs_bror->get4eBrRatio(massHiggs) / higgs_bror->getemutauBrRatio(massHiggs) * 9.;
	    	//BR_correction_2l2l = higgs_bror->get2e2muBrRatio(massHiggs) / higgs_bror->getemutauBrRatio(massHiggs) * 9. / 2.;
			BR_correction_4l = higgs_bror->get4eBrRatio(massHiggs)  * 9.;
	    	BR_correction_2l2l = higgs_bror->get2e2muBrRatio(massHiggs)  * 9. / 2.;
		}
		// Choosing the branching ratio
		if(sampleProdType == sampleType::Background || sampleProdType == sampleType::ggF_ZpZp) BR_correction = 1;
		else if(higgs->getQuadTypeBR() == quadType::Mu4 || higgs->getQuadTypeBR() == quadType::El4) BR_correction = BR_correction_4l;
		else if(higgs->getQuadTypeBR() == quadType::Mu2El2 || higgs->getQuadTypeBR() == quadType::El2Mu2) BR_correction = BR_correction_2l2l;

		// For inclusive Samples, br = 1

	
		// Dividing Lumi to match the units on the cross section...
		lumi = pileupTool->GetIntegratedLumi()/1000;
	
	}
	// Nothing for data
	if(!isMC)
	{
		crossSection = 1;
		BR_correction = 1;
		lumi = 1;
	}

	higgs->crossSection = crossSection;
	higgs->branchRatio = BR_correction;
	higgs->lumi = lumi;

	if(isDebugCall)
	{
		cout<<"------------------------"<<endl;
		cout<<"CrossSection: "<<crossSection<<endl;
		cout<<"branchRatio: "<<BR_correction<<endl;
		cout<<"lumi: "<<lumi<<endl;
		cout<<"------------------------"<<endl;
		
	}
	//return crossSection*BR_correction*lumi;

}
// Not a perfect implementation
// Just to get lumi*xsec for normalizing plot histrogram
Double_t HiggsAnalysis::getCrossSectionWeight()
{
	// Setting up the vars
	Double_t crossSection = -1;
	Double_t BR_correction = -1;
	Double_t lumi = -1;

	if(isMC)
	{
		Int_t RunNumber = event->eventinfo.mc_channel_number();
		Double_t massHiggs = getMCHiggsMass();
		CrossSections::LHCEnergy energyLHC;
		if(dataYear == 2011) energyLHC = CrossSections::SevenTeV;
		else if(dataYear == 2012) energyLHC = CrossSections::EightTeV;
	
		// get cross section 
		if(sampleProdType == sampleType::ggF) 		crossSection = Higgs_xs.higgsprodxsecGGF(massHiggs,energyLHC);
		else if(sampleProdType == sampleType::VBF) 	crossSection = Higgs_xs.higgsprodxsecVBF(massHiggs,energyLHC);
		else if(sampleProdType == sampleType::WH)  	crossSection = Higgs_xs.higgsprodxsecWH(massHiggs,energyLHC);
		else if(sampleProdType == sampleType::ZH)  	crossSection = Higgs_xs.higgsprodxsecZH(massHiggs,energyLHC);
		else if(sampleProdType == sampleType::ttH) 	crossSection = Higgs_xs.higgsprodxsecttH(massHiggs,energyLHC);
		else if(sampleProdType == sampleType::qqF) 	crossSection = Higgs_xs.higgsprodxsecGGF(massHiggs,energyLHC);
		else if(sampleProdType == sampleType::Background) 
		{
			if(energyLHC == CrossSections::SevenTeV) crossSection = CrossSections::GetBkgCrossSection7TeV(RunNumber);
			else if (energyLHC == CrossSections::EightTeV) crossSection = CrossSections::GetBkgCrossSection8TeV(RunNumber);
		}
		
		BR_correction = higgs_bror->get4eBrRatio(massHiggs)  * 9.;
		// Choosing the branching ratio
		if(sampleProdType == sampleType::Background) BR_correction = 1;
	
		// Dividing Lumi to match the units on the cross section...
		lumi = pileupTool->GetIntegratedLumi()/1000;
	}
	// Nothing for data
	if(!isMC)
	{
		crossSection = 1;
		BR_correction = 1;
		lumi = 1;
	}

	return crossSection*BR_correction*lumi;
}


// For getting the productions catergory
void HiggsAnalysis::fillProductionChannel(QuadLepton * higgs, vector<ChargedLepton *> muObject,
		vector<ChargedLepton *> elObject, vector<ChargedLepton *> jetObjects)
{
	Bool_t nJets2 = false;
	Bool_t Mjj40_130 = false;
	Bool_t Mjj130 = false;
	Bool_t isolatedLep = false;
	Bool_t hadBDTCut = false;
	Bool_t filledLeadingJet = false;

	ChargedLepton *leadingJet = 0;
	ChargedLepton *subLeadingJet = 0;
	ChargedLepton *thirdJet = 0;
	TLorentzVector sumJets;

	// For Jet Information
	if(jetObjects.size() >= 2)
	{
		// More than 2 jets
		nJets2 = true;
		Double_t leadingJetPt = -9999;
		Double_t subLeadingJetPt = -9999;
		Double_t thirdJetPt = -9999;
		// For the leading jet
		for (Int_t i = 0; i < (Int_t) jetObjects.size(); i++)
		{
			if(jetObjects[i]->get4Momentum()->Pt() > leadingJetPt && jetObjects[i]->get4Momentum()->Pt() > 25 * 1000)
			{
				// save the leading jet
				leadingJet = jetObjects[i];
				leadingJetPt = jetObjects[i]->get4Momentum()->Pt();
				filledLeadingJet = true;
			}
		}

		// Protection against the fact that there might be no leading jet
		if(!filledLeadingJet) nJets2 = false;
		//Fill subleading only if leading is filled
		if(filledLeadingJet)
		{
			for (Int_t i = 0; i < (Int_t) jetObjects.size(); i++)
			{
				if(jetObjects[i] == leadingJet) continue;
				if(jetObjects[i]->get4Momentum()->Pt() > subLeadingJetPt )
				{
					// save the leading jet
					subLeadingJet = jetObjects[i];
					subLeadingJetPt = jetObjects[i]->get4Momentum()->Pt();
				}
			}
			sumJets = leadingJet->get4MomentumNoP() + subLeadingJet->get4MomentumNoP();

			if(sumJets.M() > 130*1000) Mjj130 = true;
			if(sumJets.M() > 40*1000 && sumJets.M() < 130*1000) Mjj40_130 = true;
		}
		if(subLeadingJet)
		{
			for (Int_t i = 0; i < (Int_t) jetObjects.size(); i++)
			{
				if(jetObjects[i] == leadingJet || jetObjects[i] == subLeadingJet) continue;
				if(jetObjects[i]->get4Momentum()->Pt() > thirdJetPt )
				{
					// save the leading jet
					thirdJet = jetObjects[i];
					thirdJetPt = jetObjects[i]->get4Momentum()->Pt();
				}
			}
		}
			
		//if(leadingJet) cout<<"JetSize: "<<jetObjects.size()<<" Leading Jet: "<<leadingJet->get4Momentum()->Pt()<<" ";
		//if(subLeadingJet) cout<<"subLeadingJet: "<<subLeadingJet->get4Momentum()->Pt()<<" ";
		//if(thridJet) cout<<"thridJet: "<<thridJet->get4Momentum()->Pt()<<" ";
		//if(leadingJet || subLeadingJet || thridJet) cout<<endl<<endl;
		//for(Int_t i = 0; i < (Int_t) jetObjects.size(); i++)
		//{
		//	// Check leading jet
		//	if(jetObjects[i]->get4Momentum()->Pt() > leadingJetPt)
		//	{
		//		// Make leading now subleading
		//		subLeadingJetPt = leadingJetPt;
		//		subLeadingJet = leadingJet;
		//		// save the leading jet
		//		leadingJet = jetObjects[i];
		//		leadingJetPt = jetObjects[i]->get4Momentum()->Pt();
		//		filledLeadingJet = true;
		//	}
		//	// Check for subleading jet
		//	else if(jetObjects[i]->get4Momentum()->Pt() > subLeadingJetPt)
		//	{
		//		// save the sub-leading jet
		//		subLeadingJet = jetObjects[i];
		//		subLeadingJetPt = jetObjects[i]->get4Momentum()->Pt();
		//	}
		//}
		//sumJets = leadingJet->get4MomentumNoP() + subLeadingJet->get4MomentumNoP();

		//if(sumJets.M() > 130*1000) Mjj130 = true;
		//if(sumJets.M() > 40*1000 && sumJets.M() < 130*1000) Mjj40_130 = true;
	}


	// For BDT output
	dijet_invmass = -999*1000;
	dijet_deltaeta = -999;
	leading_jet_pt = -999*1000;
	leading_jet_eta = -999;
	subleading_jet_pt = -999*1000;

	BDT_discriminant_VBF = -999;
	BDT_discriminant_HadVH = -999;

	// fill in the leading jets
	for(Int_t i = 0; i < (Int_t) jetObjects.size(); i++)
	{
		// Check leading jet
		//if(jetObjects[i]->get4Momentum()->Pt() > leading_jet_pt)
		if(jetObjects[i]->get4Momentum()->Pt() > leading_jet_pt && jetObjects[i]->get4Momentum()->Pt() > 25 * 1000)
		{
			leading_jet_pt = jetObjects[i]->get4Momentum()->Pt();
			leading_jet_eta = jetObjects[i]->get4Momentum()->Eta();
			leadingJet = jetObjects[i];

		}
	}
	
	if(nJets2)
	{
	
		dijet_invmass = sumJets.M();
		dijet_deltaeta = fabs(leadingJet->get4Momentum()->Eta() - subLeadingJet->get4Momentum()->Eta());
		leading_jet_pt = leadingJet->get4Momentum()->Pt();
		leading_jet_eta = leadingJet->get4Momentum()->Eta();
		subleading_jet_pt = subLeadingJet->get4Momentum()->Pt();
	
		// computing HadVH and VBF BDT discriminant and filling the variables
	   
	    BDT_discriminant_VBF   = CategoriesDiscriminantTool->Get_VBFDiscriminant_Output();
	    BDT_discriminant_HadVH = CategoriesDiscriminantTool->Get_HadVHDiscriminant_Output(); 
	
		if(CategoriesDiscriminantTool->Pass_HadVHDiscriminant() ) hadBDTCut = true;

		if(isDebugCall)
		{
			cout<<"-------------------"<<endl;
			cout<<"dijet_invmass: "<<dijet_invmass<<endl;
			cout<<"dijet_deltaeta: "<<dijet_deltaeta<<endl;
			cout<<"leading_jet_pt: "<<leading_jet_pt<<endl;
			cout<<"leading_jet_eta: "<<leading_jet_eta<<endl;
			cout<<"subleading_jet_pt: "<<subleading_jet_pt<<endl;
			cout<<"dijet_invmass: "<<dijet_invmass<<endl;
			cout<<"BDT_discriminant_VBF: "<<BDT_discriminant_VBF<<endl;
			cout<<"BDT_discriminant_HadVH: "<<BDT_discriminant_HadVH<<endl;
			cout<<"-------------------"<<endl;			
		}


	}
	vector<ChargedLepton *> extraLep;
	// For isolated Lepton
	if((muObject.size() + elObject.size()) > 4)
	{
		for(Int_t i = 0; i < (Int_t) muObject.size(); i++)
		{
			if(isGoodExtraLepton(higgs, muObject[i]))
			{
				//higgs->SetProductionChannel(productionChannel::VH);
				//return;
				isolatedLep = true;
				fillExtraLepParent(muObject[i]);
				extraLep.push_back(muObject[i]);
			}
		}
		for(Int_t i = 0; i < (Int_t) elObject.size(); i++)
		{
			if(isGoodExtraLepton(higgs, elObject[i]))
			{
				//higgs->SetProductionChannel(productionChannel::VH);
				//return;
				isolatedLep = true;
				fillExtraLepParent(elObject[i]);				
				extraLep.push_back(elObject[i]);				
			}
		}
	}
	if(isDebugCall) cout<<"Size of extraLep: "<<extraLep.size()<<endl;
	
	ChargedLepton *leadingExtraLep = 0;
	ChargedLepton *subleadingExtraLep = 0;
	Double_t lep_leadingPt = -999999;
	Double_t sublep_leadingPt = -999999;
	for(Int_t i = 0; i < (Int_t) extraLep.size(); i++)
	{
		if(extraLep[i]->get4Momentum()->Pt() > lep_leadingPt)
		{
			sublep_leadingPt = lep_leadingPt;
			subleadingExtraLep = leadingExtraLep;
			
			leadingExtraLep = extraLep[i];
			lep_leadingPt = extraLep[i]->get4Momentum()->Pt();
		}
		else if(extraLep[i]->get4Momentum()->Pt() > sublep_leadingPt)
		{
			subleadingExtraLep = extraLep[i];
			sublep_leadingPt = extraLep[i]->get4Momentum()->Pt();
		}
	}

	// For ggF catergorization
	//higgs->SetProductionChannel(productionChannel::ggF);

	//Setting the categories
	if(nJets2 && Mjj130) higgs->SetProductionChannel(productionChannel::VBF);
	else if(nJets2 && Mjj40_130 && hadBDTCut) higgs->SetProductionChannel(productionChannel::VHHad);
	else if(isolatedLep) higgs->SetProductionChannel(productionChannel::VHLep);
	else higgs->SetProductionChannel(productionChannel::ggF);

	//Filling variables for the tree
	higgs->leadingJet 				= leadingJet;
	higgs->subLeadingJet 			= subLeadingJet;
	higgs->thirdJet			 		= thirdJet;
	higgs->n_jets 					= jetObjects.size();
	higgs->dijet_invmass		 	= dijet_invmass/1000;
	higgs->dijet_deltaeta 			= dijet_deltaeta;
	higgs->leading_jet_pt 			= leading_jet_pt/1000;
	higgs->leading_jet_eta 			= leading_jet_eta;
	higgs->subleading_jet_pt 		= subleading_jet_pt/1000;
	higgs->BDT_discriminant_VBF		= BDT_discriminant_VBF;
	higgs->BDT_discriminant_HadVH 	= BDT_discriminant_HadVH;
	higgs->leadingExtraLep			= leadingExtraLep;
	higgs->subleadingExtraLep		= subleadingExtraLep;

	if(isDebugCall)
	{
		cout<<"-----------"<<endl<<"Lepton type info"<<endl;
		if(leadingExtraLep != 0)
		{
			cout<<"leading letpton"<<endl;
			cout<<"pT: "<<leadingExtraLep->get4Momentum()->Pt()<<endl;
			cout<<"M: "<<leadingExtraLep->get4Momentum()->M()<<endl;
			cout<<"LepType: "<<leadingExtraLep->lepType<<endl;
			cout<<"Parent LepType: "<<leadingExtraLep->truthParentType<<endl;
			cout<<"Charge: "<<leadingExtraLep->charge<<endl;
		}
		if(subleadingExtraLep != 0)
		{
			cout<<"subleadingExtraLep letpton"<<endl;
			cout<<"pT: "<<subleadingExtraLep->get4Momentum()->Pt()<<endl;
			cout<<"M: "<<subleadingExtraLep->get4Momentum()->M()<<endl;
			cout<<"LepType: "<<subleadingExtraLep->lepType<<endl;
			cout<<"Parent LepType: "<<subleadingExtraLep->truthParentType<<endl;
			cout<<"Charge: "<<subleadingExtraLep->charge<<endl;
			
		}

		cout<<endl;
	}
	extraLep.clear();
	
}
// To fill in for fudicial jets
void HiggsAnalysis::fillFudicialJets(QuadLepton * higgs)
{
	Double_t leading_jet_pt = -999*1000;
	// fill in the leading jets
	for(Int_t i = 0; i < (Int_t) jetsEvent_Fid.size(); i++)
	{
		// Check leading jet
		if(jetsEvent_Fid[i]->get4Momentum()->Pt() > leading_jet_pt)
		{
			leading_jet_pt = jetsEvent_Fid[i]->get4Momentum()->Pt();

		}
	}
	Double_t leading_jet_pt_truth = -999*1000;
	// fill in the leading jets
	for(Int_t i = 0; i < (Int_t) jetsTruthEvent_Fid.size(); i++)
	{
		// Check leading jet
		if(jetsTruthEvent_Fid[i]->get4Momentum()->Pt() > leading_jet_pt_truth)
		{
			leading_jet_pt_truth = jetsTruthEvent_Fid[i]->get4Momentum()->Pt();

		}
	}

	higgs->n_jets_fid = jetsEvent_Fid.size();
	higgs->leading_jet_pt_fid = leading_jet_pt/1000;
	higgs->n_jets_truth_fid = jetsTruthEvent_Fid.size();
	higgs->leading_jet_pt_truth_fid = leading_jet_pt_truth/1000;

	//cout<<"Njets_fid: "<<jetsEvent_Fid.size()<<" pT: "<<leading_jet_pt<<endl;
	//cout<<"Njets_fid_Truth: "<<jetsTruthEvent_Fid.size()<<" pT: "<<leading_jet_pt_truth<<endl;

}

// To find addition extra good lepton for VH catergorization
Bool_t HiggsAnalysis::isGoodExtraLepton(QuadLepton * higgs, ChargedLepton * lep)
{
	if(isDebugCall && !isDebugCallllee) cout<<"---------------"<<endl<<"checking extra lepton"<<endl;	
	vector <ChargedLepton *> higgsLep = higgs->getLepton();
	// Check for overlap
	for(vector<ChargedLepton *>::iterator itr_lep = higgsLep.begin();
		itr_lep != higgsLep.end(); ++itr_lep)
		{
			ChargedLepton* curr_i = *itr_lep;
			// if same return false
			if(curr_i == lep) return false;
		}

	// Kinematic Cut
	// return false if pT < 8Gev
	if(isDebugCall) cout<<"pass pT cut"<<endl;
	if(!(lep->get4Momentum()->Pt() > 8*1000)) return false;

	Bool_t passCut = true;
	// deltaR cut
	for(vector<ChargedLepton *>::iterator itr_i = higgsLep.begin();
		itr_i != higgsLep.end(); ++itr_i)
	{
					// Getting the charged Lepton
			ChargedLepton * lep_i = *itr_i;
			ChargedLepton * lep_j = lep;
		
			// Since the cut value is dependant on flavour type
			// getting that
			Double_t cutDeltaR = 0;
			if(lep_i->getFlavor() == lep_j->getFlavor()) {cutDeltaR = 0.10;}
			else {cutDeltaR = 0.20;}

			// Actually compute the delta R and cut on it
			Double_t dataDeltaR = 0;
			dataDeltaR = lep_i->get4Momentum()->DeltaR(*(lep_j->get4Momentum()));
			
			// Debugging information
			if(isDebugCall){cout<<"Delta R "<<dataDeltaR<<"\t";}

			if(dataDeltaR < cutDeltaR) passCut = false;

			//if(!passCut) break;
			// Debugging information
			if(isDebugCall){cout<<endl;}
	}
	if(!passCut) return false;
	if(isDebugCall) cout<<"EL: pass deltaR"<<endl;

	// Track iso
	Double_t cutDeltaR = 0.20;
	Double_t curr_PTCone20 = lep->getPTCone20();
	
	for(vector<ChargedLepton *>::iterator itr_j = higgsLep.begin(); itr_j != higgsLep.end(); itr_j++)
	{
		ChargedLepton *lep_j = *itr_j;

		Double_t dataDeltaR = lep->get4Momentum()->DeltaR(*(lep_j->get4Momentum()));
		Double_t toSubtract = 0;
		if(dataDeltaR < cutDeltaR) 
		{
			if(lep_j->getType() == leptonType::MuonStaco || lep_j->getType() == leptonType::MuonCalo)
			{
				toSubtract = (1/fabs((lep_j->GetMuon()->id_qoverp()))*sin((lep_j->GetMuon()->id_theta())));
			}
			else if (lep_j->getType() == leptonType::MuonStandAlone)
			{
				toSubtract = 0;
			}
			else if(lep_j->getType() == leptonType::ElectronGSF)
			{
				toSubtract = lep_j->GetElectron()->trackpt();
			}
		}
		curr_PTCone20 = curr_PTCone20 - toSubtract;
	}
	
	Double_t dataTrackIso = curr_PTCone20/lep->get4Momentum()->Pt();
	if(dataTrackIso > 0.15) return false;
	if(isDebugCall) cout<<"pass track iso cut"<<endl;

	// Calo Iso
	Double_t EtCone20CorrectedIso = 0;
	if(lep->getFlavor() == flavor::Electron)
	{
		D3PDReader::ElectronD3PDObjectElement* el_curr = lep->GetElectron();
		if(dataYear == 2011)
		{
			EtCone20CorrectedIso = el_curr->Etcone20();
		}
		else if(dataYear == 2012)
		{
			EtCone20CorrectedIso = el_curr->topoEtcone20();
		}
		// Smearing Correction
		if(doCaloIsolationCorr && doCorr) EtCone20CorrectedIso = getCaloIsoCorrection(el_curr);
	}
	else if(lep->getFlavor() == flavor::Muon)
	{
		D3PDReader::MuonD3PDObjectElement* mu_curr = lep->GetMuon();
		EtCone20CorrectedIso = mu_curr->etcone20();
	}
		if(isDebugCall) cout<<"EL: OriginalIso "<<EtCone20CorrectedIso<<endl;

	// This gets what I need to subtract in the end for calo iso
	for(vector<ChargedLepton *>::iterator itr_i = higgsLep.begin(); itr_i != higgsLep.end(); itr_i++)
	{
		// Just to keep track of what event loop we are at
		ChargedLepton *lep_i = *itr_i;
		Double_t EtOverlap = 0;
		if(lep_i->getFlavor() == flavor::Electron)
		{
			D3PDReader::ElectronD3PDObjectElement* el_curr = lep_i->GetElectron();
			EtOverlap = el_curr->rawcl_E()/cosh(el_curr->tracketa());
		}
		else if(lep_i->getFlavor() == flavor::Muon)
		{
			D3PDReader::MuonD3PDObjectElement* mu_curr = lep_i->GetMuon();
			EtOverlap = mu_curr->pt();
			EtOverlap = 0;
		}
		if(isDebugCall) cout<<"EL: ToSubtract: EtOverlap "<<EtOverlap<<endl;
		if(DeltaRCaloIso(lep, lep_i) < 0.18)
		{
			EtCone20CorrectedIso = EtCone20CorrectedIso - EtOverlap;
			if(isDebugCall) cout<<"EL: passdR sub:"<<EtCone20CorrectedIso<<endl;
		}
	}
	// Cut Values;
	Double_t muonCut = 0.30;
	Double_t muonStandAloneCut = 0.15;

	Double_t electronCut2011 = 0.30;
	Double_t electronCut2012 = 0.20;

	Double_t dataCut = EtCone20CorrectedIso/lep->get4Momentum()->Pt();
	if(isDebugCall)
	{
		cout<<"EL: EtCone20CorrectedIso: "<<EtCone20CorrectedIso<<endl;
		cout<<"EL: pT: "<<lep->get4Momentum()->Pt()<<endl;
		cout<<"EL: dataCut: "<<dataCut<<endl;
	}
	if(lep->getFlavor() == flavor::Electron)
	{
		if(dataYear == 2011) 
		{
			if(dataCut > electronCut2011) return false;
		}
		else if(dataYear == 2012)
		{
			if(dataCut > electronCut2012) return false;
		}
	}
	if(lep->getFlavor() == flavor::Muon)
	{
		if(lep->getType() != leptonType::MuonStandAlone) 
		{
			if(dataCut > muonCut) return false;
		}
		else if(lep->getType() == leptonType::MuonStandAlone)
		{
			if(dataCut > muonStandAloneCut) return false;
		}
	}
	if(isDebugCall) cout<<"pass calo iso cut"<<endl;

	// DO Sig
	// Cut Values;
	Double_t muD0SigCut = 3.5;
	Double_t elD0SigCut = 6.5;
	
	// Actual Cuts
	if(lep->getFlavor() == flavor::Electron)
	{
		Double_t dataD0Sig = fabs(lep->GetElectron()->trackd0pvunbiased())/ lep->GetElectron()->tracksigd0pvunbiased();
		if(dataD0Sig > elD0SigCut) return false;
	}
	if(lep->getFlavor() == flavor::Muon)
	{
		Double_t dataD0Sig = fabs(lep->GetMuon()->trackd0pvunbiased())/ lep->GetMuon()->tracksigd0pvunbiased();
		if(dataD0Sig > muD0SigCut) return false;			
	}
	if(isDebugCall) cout<<"pass d0 iso cut"<<endl;
	

	// Calo/Standalone count
	Int_t caloSACount = 0;
	if(lep->getType() == leptonType::MuonCalo || lep->getType() == leptonType::MuonStandAlone) caloSACount++;
	for(vector<ChargedLepton *>::iterator itr_i = higgsLep.begin(); itr_i != higgsLep.end(); itr_i++)
	{
		// Just to keep track of what event loop we are at
		ChargedLepton *lep_i = *itr_i;
		if(lep_i->getType() == leptonType::MuonCalo || lep_i->getType() == leptonType::MuonStandAlone) caloSACount++;
	}

	if(caloSACount > 1) return false;
	if(isDebugCall) cout<<"pass calo and standAlone cut"<<endl<<"----------------"<<endl;


	return true;
}
void HiggsAnalysis::fillExtraLepParent(ChargedLepton* lep)
{
	Int_t lepBareTruthIndex = -1;
	if(lep->getFlavor() == flavor::Electron)
	{
		D3PDReader::ElectronD3PDObjectElement* el =  lep->GetElectron();
		// Get the barCode
		Int_t truthBarcode = el->truth_barcode();
		Int_t truthPDG = el->truth_type();
	
		if(abs(truthPDG) != 11) {truthBarcode = -1;}
		
		// Pushing in the truth barcode
		lepBareTruthIndex = getIndexBarcodeMatch(truthBarcode, truthPDG);
	}
	else if(lep->getFlavor() == flavor::Muon)
	{
		D3PDReader::MuonD3PDObjectElement* mu =  lep->GetMuon();
		Int_t truthBarcode = mu->truth_barcode();
		Int_t truthPDG = mu->truth_type();
		
		lepBareTruthIndex = getIndexBarcodeMatch(truthBarcode, truthPDG);
	}

	// Get the Z index and mass
	vector<Int_t> zHiggsIndex;
	vector<Double_t> zHiggsMass;
	for (Int_t i = 0; i < (Int_t) event->mc.n(); i++)
	{
		Int_t pdgid = event->mc[i].pdgId();
		Int_t status = event->mc[i].status();
		
		// Look for a Z that decays
		if(pdgid == 23 &&  ((status==22    || // Z->4l (Z coming from H->ZZ                     (l status 3 and 1)
                            status==155   || // Jimmy gg2ZZ NOT WORKING
                            status==10902 || // Pythia6 stable Z->ll (equivalent to status==2) (l status 1)
                            status==2     || // Pythia6 stable Z->ll                           (l status 1)
                            status==3     || // Pythia6 stable Z->ll                           (l status 3)
                            status==62)))    // PowHegPythia8 ZZ->4l                           (l status 3 and 1))
		{
			// To veto the third Z in ZH samples
			Bool_t higgsDecayZ = true;
			if(sampleProdType == sampleType::ZH) higgsDecayZ = checkParentHiggs(i);
			if(!higgsDecayZ) continue;

			zHiggsIndex.push_back(i);
			zHiggsMass.push_back(event->mc[i].m());

		}
	}

	if(zHiggsIndex.size() != 2 && generatorName == MCGeneratorName::Pythia && sampleProdType != sampleType::Background)
		cout<<"Error: too many (or too less) Z boson for extra lep: size: "<<zHiggsIndex.size()<<endl;

	Int_t z1Index = -1;
	Int_t z2Index = -1;

	if(zHiggsIndex.size() == 2)
	{
		if(zHiggsMass[0] > zHiggsMass[1])
		{
			z1Index = zHiggsIndex[0];
			z2Index = zHiggsIndex[1];
		}
		else
		{
			z1Index = zHiggsIndex[1];
			z2Index = zHiggsIndex[0];
		}

	}

	// So have the truth bare index
	// Now check if the parents are W
	Bool_t parentW = checkParent(lepBareTruthIndex, 24) || checkParent(lepBareTruthIndex, -24);

	Bool_t parentZ = checkParent(lepBareTruthIndex, 23);
	Bool_t parentZ1 = checkParent(lepBareTruthIndex, 23, z1Index);
	Bool_t parentZ2 = checkParent(lepBareTruthIndex, 23, z2Index);
	Bool_t parentHiggs = checkParent(lepBareTruthIndex, 25) || checkParent(lepBareTruthIndex, 39);

	zHiggsIndex.clear();
	zHiggsMass.clear();

	if(parentW) lep->truthParentType = VHTruthType::WnoHiggs;
	else if(parentZ && !parentHiggs) lep->truthParentType = VHTruthType::ZnoHiggs;
	else if(parentZ && parentHiggs && parentZ1) lep->truthParentType = VHTruthType::Z1fromHiggs;
	else if(parentZ && parentHiggs && parentZ2) lep->truthParentType = VHTruthType::Z2fromHiggs;
	else lep->truthParentType = VHTruthType::unknown;



}


// calculates and fills the BDT weight
void HiggsAnalysis::fillBDTWeights(QuadLepton * higgs)
{

	// For pt Systematics
	Float_t pt4l_truth_born = -999;
	if(isMC && higgs->truthVec.Pt() != 0) pt4l_truth_born = higgs->truthVec.Pt()/1000; 

	Bool_t isggH12 = false;
	if(isMC && dataYear == 2012 && sampleProdType == sampleType::ggF) isggH12 = true;

	// Type of quad we are dealing with...
	H4lBDTWeights::BDTQuadType BDTquad = H4lBDTWeights::UNKNOWN;

	if(higgs->getQuadType() == quadType::Mu4) BDTquad = H4lBDTWeights::_4mu;
	else if (higgs->getQuadType() == quadType::Mu2El2) BDTquad = H4lBDTWeights::_2mu2e;
	else if (higgs->getQuadType() == quadType::El2Mu2) BDTquad = H4lBDTWeights::_2e2mu;
	else if(higgs->getQuadType() == quadType::El4) BDTquad = H4lBDTWeights::_4e;


	// Give Input to the BDT tool
	BDTtool->setBDTInputs(higgs->getZ1()->getLepPlus()->get4MomentumBDT(),
			   		   	  higgs->getZ1()->getLepNeg()->get4MomentumBDT(),
					      higgs->getZ2()->getLepPlus()->get4MomentumBDT(),
			   		      higgs->getZ2()->getLepNeg()->get4MomentumBDT(),
					      higgs->sum_unconstrained.Pt()/1000,
					      higgs->sum_unconstrained.Eta(),
					      pt4l_truth_born,
					      BDTquad,
					      isggH12);

	// Output	
	BDTtool->fillBDTOutputs(higgs->KD_discriminant, 
							higgs->BDT_discriminant,
							higgs->BDTGuass_discriminant, 
							higgs->ptSysupFac, 
							higgs->ptSysdownFac);
	//cout.setf(ios::fixed, ios::floatfield);
	//cout.setf(ios::showpoint);

	if(isDebugCall)
	{
		cout<<"-------------------------"<<endl;
		TLorentzVector temp = higgs->getZ1()->getLepPlus()->get4MomentumBDT();
		cout<<"Z1Lep Plus Pt: "<<temp.Pt()<<" eta: "<<temp.Eta()<<" Phi: "<<temp.Phi()<<" M: "<<temp.M()<<" E: "<<temp.E()<<endl;
		temp = higgs->getZ1()->getLepNeg()->get4MomentumBDT();
		cout<<"Z1Lep Neg Pt: "<<temp.Pt()<<" eta: "<<temp.Eta()<<" Phi: "<<temp.Phi()<<" M: "<<temp.M()<<" E: "<<temp.E()<<endl;
		temp = higgs->getZ2()->getLepPlus()->get4MomentumBDT();
		cout<<"Z2Lep Plus Pt: "<<temp.Pt()<<" eta: "<<temp.Eta()<<" Phi: "<<temp.Phi()<<" M: "<<temp.M()<<" E: "<<temp.E()<<endl;
		temp = higgs->getZ2()->getLepNeg()->get4MomentumBDT();
		cout<<"Z2Lep Neg Pt: "<<temp.Pt()<<" eta: "<<temp.Eta()<<" Phi: "<<temp.Phi()<<" M: "<<temp.M()<<" E: "<<temp.E()<<endl;

		cout<<"pT_unconstrained: "<<higgs->sum_unconstrained.Pt()/1000<<endl;
		cout<<"Eta_unconstrained: "<<higgs->sum_unconstrained.Eta()<<endl;
		cout<<"pt4l_truth_born: "<<pt4l_truth_born<<endl;
		cout<<"BDTquad: "<<BDTquad<<endl;
		cout<<"isggH12: "<<isggH12<<endl<<endl;

		cout<<"KD_discriminant: "<<higgs->KD_discriminant<<endl;
		cout<<"BDT_discriminant: "<<higgs->BDT_discriminant<<endl;
		cout<<"BDTGuass_discriminant: "<<higgs->BDTGuass_discriminant<<endl;
		cout<<"ptSysupFac: "<<higgs->ptSysupFac<<endl;
		cout<<"ptSysdownFac: "<<higgs->ptSysdownFac<<endl;
		cout<<"-------------------------"<<endl;
	}
}

// Decides on the analysis based on the data stream name 
void HiggsAnalysis::fillStreamAnalysis()
{
	TString currentFile;
	currentFile = currFileName;
	// Overwrite from grid
	if(runningGrid) currentFile = gridFileName;

	if(currentFile.Contains("_Egamma"))
	{
		streamName = streamContainer::eGamma;
		do4Mu = false;
		do4El = true;
		do2L2L = true;
	}
	else if(currentFile.Contains("_Muons"))
	{
		streamName = streamContainer::Muon;
		do4Mu = true;
		do4El = false;
		do2L2L = true;
	}
	else
	{
		streamName = streamContainer::Other;
		do4Mu = true;
		do4El = true;
		do2L2L = true;
	}
}
// Flags the events based on comparing the higest pT lepton and stream
void HiggsAnalysis::fillFlagStream(QuadLepton * higgs)
{
	vector <ChargedLepton *> higgsLep = higgs->getLepton();
	
	Double_t maxPt = -99999;
	Int_t maxPtFlavor= -1;

	for(vector<ChargedLepton *>::iterator itr_lep = higgsLep.begin();
		itr_lep != higgsLep.end(); ++itr_lep)
		{
			ChargedLepton* curr_i = *itr_lep;
			// find the highest pT lepton
			if(curr_i->get4Momentum()->Pt() > maxPt)
			{
				maxPt = curr_i->get4Momentum()->Pt();
				maxPtFlavor = curr_i->getFlavor();
			}
		}
	// Sanity Check
	if(maxPt < 0) cout<<"ERROR: fillFlagStream, higestPT algorithm not working"<<endl;
	// Checking for flagging
	Bool_t flagQuad = false;
	if(streamName == streamContainer::Muon && maxPtFlavor != flavor::Muon) flagQuad = true;
	else if(streamName == streamContainer::eGamma && maxPtFlavor != flavor::Electron) flagQuad = true;
	else flagQuad = false;

	// Saving the result
	higgs->flagQuad = (Int_t) flagQuad;	

}

// Fill truth information for any lepton
void HiggsAnalysis::fillTruthRecoMatchedInfo(QuadLepton * higgs)
{
	// To safeguard against data...
	if(!isMC) return;
	// Getting the information on barcode and the subsequent
	vector<ChargedLepton *> lep_i = higgs->getLepton();
	vector<Int_t> lepBareTruthIndex;
	vector<Int_t> lepBareTruthPDGID;
	vector<TLorentzVector> lep_bareVec;
	
	Int_t j = -1;
	for(vector<ChargedLepton *>::iterator itr_lep = lep_i.begin();
		itr_lep != lep_i.end(); ++itr_lep)
	{
		// Initializing the variables
		ChargedLepton* curr_i = *itr_lep;
		lepBareTruthIndex.push_back(-1);
		lepBareTruthPDGID.push_back(-1);

		Int_t truthBarcode = -1;
		Int_t truthPDG = -1;
		j++;

		if(curr_i->getFlavor() == flavor::Electron)
		{
			D3PDReader::ElectronD3PDObjectElement* el =  curr_i->GetElectron();
			// Get the barCode
			truthBarcode = el->truth_barcode();
			truthPDG = el->truth_type();

			if(abs(truthPDG) != 11) {truthBarcode = -1;}

			// Ensure that the mother is a Z, otherwise this is a fake. Specific to Pythia
			if(abs(el->truth_mothertype()) != 23 && abs(el->truth_mothertype()) != 15) {truthBarcode = -1;}

			// Pushing in the truth barcode
			lepBareTruthIndex[j] = getIndexBarcodeMatch(truthBarcode, truthPDG);
			lepBareTruthPDGID[j] = truthPDG;
		}
		else if(curr_i->getFlavor() == flavor::Muon)
		{
			D3PDReader::MuonD3PDObjectElement* mu =  curr_i->GetMuon();
			truthBarcode = mu->truth_barcode();
			truthPDG = mu->truth_type();

			// Fake
			if(abs(truthPDG) != 13) {truthBarcode = -1;}
			// Ensure that the mother is a Z, otherwise this is a fake. Specific to Pythia
			if(abs(mu->truth_mothertype()) != 23 && abs(mu->truth_mothertype()) != 15) {truthBarcode = -1;}
			lepBareTruthIndex[j] = getIndexBarcodeMatch(truthBarcode, truthPDG);

			if(curr_i->getType() == leptonType::MuonStandAlone) 
			{
				if(curr_i->charge > 0) truthPDG = 13;
				else truthPDG = -13;

				// Current hack
				truthPDG = fabs(truthPDG);
				lepBareTruthIndex[j] = getIndexDeltaRMatch(truthPDG, curr_i->get4Momentum()->Eta(), curr_i->get4Momentum()->Phi());
				if(lepBareTruthIndex[j] > 0) truthPDG = event->mc[lepBareTruthIndex[j]].pdgId();
			}
			lepBareTruthPDGID[j] = truthPDG;
		}
	}
	if(isDebugCall) cout<<"---------------------"<<endl<<"Parent info for Lepton"<<endl;
	// Checking if the parents of the leptons are Higgs and Z
	for(Int_t i = 0; i < (Int_t) lepBareTruthIndex.size(); i++)
	{
		Bool_t parentZ = checkParent(lepBareTruthIndex[i], 23);
		Bool_t parentHiggs = checkParent(lepBareTruthIndex[i], 25) || checkParent(lepBareTruthIndex[i], 39);
		if(isDebugCall) 
		{
			cout<<"Lep Index: "<<lepBareTruthIndex[i]<<" has Z parent: ";
			cout<<parentZ;
			cout<<" has parent Higgs: ";
			cout<<parentHiggs<<endl;
		}

		// to get reco info for bkg samples.
		if(sampleProdType == sampleType::Background) parentHiggs = true;

		if(!(parentZ && parentHiggs)) lepBareTruthIndex[i] = -1;
	}

	// Check the size of the lepton
	if(lepBareTruthIndex.size() != 4) cout<<"Error: fillTruthRecoMatchedInfo: Too many leptons"<<endl;
	lep_bareVec = getVecFromIndex(lepBareTruthIndex);

	// Check the size of the lepton
	if(lep_bareVec.size() != 4) cout<<"Error: fillTruthRecoMatchedInfo: Too many lep_bareVec"<<endl;

	// Now fill in the bare Z1, Z2 and the Higgs Mass
	// Due to the original numbering of the lepton in the higgs, 
	// 0,1 are lepplus and lepminus for z1
	// 2,3 are lepplus and lepminus for z2
	Bool_t fillZ1Bare = true;
	Bool_t fillZ2Bare = true;

	if(lep_bareVec[0].Pt() == 0) fillZ1Bare = false;
	if(lep_bareVec[1].Pt() == 0) fillZ1Bare = false;
	if(lep_bareVec[2].Pt() == 0) fillZ2Bare = false;
	if(lep_bareVec[3].Pt() == 0) fillZ2Bare = false;
	
	if(fillZ1Bare) higgs->mZ1_truth_bare = (lep_bareVec[0] + lep_bareVec[1]).M()/1000;	
	if(fillZ2Bare) higgs->mZ2_truth_bare = (lep_bareVec[2] + lep_bareVec[3]).M()/1000;
	if(fillZ1Bare && fillZ2Bare) higgs->m4l_truth_bare = (lep_bareVec[0] + lep_bareVec[1] + lep_bareVec[2] + lep_bareVec[3]).M()/1000;

	// Born
	// This step is now specific to pythia*, and powheg. Will not work for sherpa or Jimmy
	vector<Int_t> lepBornTruthIndex = getBornIndexFromBare(lepBareTruthIndex, lepBareTruthPDGID, true);
	vector<Int_t> lepMotherTruthIndex = getMotherIndexFromBare(lepBareTruthIndex, lepBareTruthPDGID, lepBornTruthIndex);
	
	// Filling the born TLorentzVector and 
	// Now that we have the status 1 leptons, we can fill the TLorentzVector for the bare 
	vector<TLorentzVector> lep_bornVec;
	lep_bornVec = getVecFromIndex(lepBornTruthIndex);
	
	// Now fill in the bare Z1, Z2 and the Higgs Mass
	// Due to the original numbering of the lepton in the higgs, 
	// 0,1 are lepplus and lepminus for z1
	// 2,3 are lepplus and lepminus for z2
	Bool_t fillZ1Born = true;
	Bool_t fillZ2Born = true;
	if(lepBornTruthIndex.size() != 4) cout<<"lepBornTruthIndex size: "<<lepBornTruthIndex.size() <<endl;
	
	if(lep_bornVec[0].Pt() == 0) fillZ1Born = false;
	if(lep_bornVec[1].Pt() == 0) fillZ1Born = false;
	if(lep_bornVec[2].Pt() == 0) fillZ2Born = false;
	if(lep_bornVec[3].Pt() == 0) fillZ2Born = false;

	if(fillZ1Born) higgs->mZ1_truth_born = (lep_bornVec[0] + lep_bornVec[1]).M()/1000;
	if(fillZ2Born) higgs->mZ2_truth_born = (lep_bornVec[2] + lep_bornVec[3]).M()/1000;
	if(fillZ1Born && fillZ2Born) higgs->m4l_truth_born = (lep_bornVec[0] + lep_bornVec[1] + lep_bornVec[2] + lep_bornVec[3]).M()/1000;

	// Cleaning the parent information
	std::sort(lepMotherTruthIndex.begin(), lepMotherTruthIndex.end());	// 
	lepMotherTruthIndex.erase(std::unique(lepMotherTruthIndex.begin(), lepMotherTruthIndex.end()), lepMotherTruthIndex.end());
	
	// Getting the FSR Photon
	// Looping over the MC to find photon with status 1 and the same Z parents
	vector<Int_t> lepFSRTruthIndex;	
	vector<TLorentzVector> lepFSRTruthVec;	
	
	for (Int_t i = 0; i < (Int_t) event->mc.n(); i++)
	{
		// Check for type and status
		if(event->mc[i].pdgId() == 22 && event->mc[i].status() == 1)
		{
			if(event->mc[i].parent_index().size() == 0 ) continue;
			// Looping over parent to check it the parent in one of the Z
			for(Int_t j = 0; j < (Int_t) event->mc[i].parent_index().size(); j++)
			{
				if(std::find(lepMotherTruthIndex.begin(), lepMotherTruthIndex.end(), event->mc[i].parent_index().at(j)) != lepMotherTruthIndex.end())
				{
					lepFSRTruthIndex.push_back(i);
					break;
				}
			}
		}
	}

	// Clean FSR photon
	std::sort(lepFSRTruthIndex.begin(), lepFSRTruthIndex.end());
	lepFSRTruthIndex.erase(std::unique(lepFSRTruthIndex.begin(), lepFSRTruthIndex.end()), lepFSRTruthIndex.end());

	// Create TLorentzVector for bare zBoson
	TLorentzVector vecZ1; 
	TLorentzVector vecZ2;
	vecZ1.SetPtEtaPhiM(0,0,0,0);
	vecZ2.SetPtEtaPhiM(0,0,0,0);
	// For Debug counting
	Int_t ndRPassZ1 = 0;
	Int_t ndRPassZ2 = 0;

	if(fillZ1Bare) vecZ1 =  lep_bareVec[0] + lep_bareVec[1];
	if(fillZ2Bare) vecZ2 =  lep_bareVec[2] + lep_bareVec[3];
	// Computing delta R for FSR and adding them to the right photon
	for(Int_t i = 0; i < (Int_t) lepFSRTruthIndex.size(); i++)
	{
		Int_t currI = lepFSRTruthIndex[i];
		TLorentzVector FSRPhoton;
		FSRPhoton.SetPtEtaPhiM(event->mc[currI].pt(), 
							  event->mc[currI].eta(), 
							  event->mc[currI].phi(), 
							  event->mc[currI].m());
		lepFSRTruthVec.push_back(FSRPhoton);

		// Loop over z1
		for(Int_t j = 0; j <= 1; j++)
		{
			if(vecZ1.Pt() == 0 || lep_bareVec[j].Pt() == 0) continue;	
			// Check if the same parent
			if(!checkSameParent(currI, lepBareTruthIndex[j])) continue;
			if(lep_bareVec[j].DeltaR(FSRPhoton) < 0.1) {vecZ1 = vecZ1 + FSRPhoton;ndRPassZ1++; break;}
		}
		// Loop over z2
		for(Int_t j = 2; j <= 3; j++)
		{
			if(vecZ2.Pt() == 0 || lep_bareVec[j].Pt() == 0) continue;
			// Check if the same parent
			if(!checkSameParent(currI, lepBareTruthIndex[j])) continue;
			if(lep_bareVec[j].DeltaR(FSRPhoton) < 0.1 ) {vecZ2 = vecZ2 + FSRPhoton;ndRPassZ2++; break;}
			
		}
	}
	
	if(fillZ1Bare) higgs->mZ1_truth_dressed = vecZ1.M()/1000;
	if(fillZ2Bare) higgs->mZ2_truth_dressed = vecZ2.M()/1000;
	if(fillZ1Bare && fillZ2Bare) higgs->m4l_truth_dressed = (vecZ1 + vecZ2).M()/1000;

	// Saving the vectors for later information
	higgs->getZ1()->getLepPlus()->m_momentumTruthRecoBare	= lep_bareVec[0];	 
	higgs->getZ1()->getLepNeg()->m_momentumTruthRecoBare	= lep_bareVec[1];
	higgs->getZ2()->getLepPlus()->m_momentumTruthRecoBare	= lep_bareVec[2];
	higgs->getZ2()->getLepNeg()->m_momentumTruthRecoBare 	= lep_bareVec[3];

	higgs->getZ1()->getLepPlus()->m_momentumTruthRecoBorn	= lep_bornVec[0];
	higgs->getZ1()->getLepNeg()->m_momentumTruthRecoBorn	= lep_bornVec[1];
	higgs->getZ2()->getLepPlus()->m_momentumTruthRecoBorn	= lep_bornVec[2];
	higgs->getZ2()->getLepNeg()->m_momentumTruthRecoBorn 	= lep_bornVec[3];


	// Filling the ture-truth information
	TLorentzVector bornSumVec = (lep_bornVec[0] + lep_bornVec[1] + lep_bornVec[2] + lep_bornVec[3]);

	if(isDebugCall)
	{
		cout<<"--------------------------------"<<endl;
		cout<<"Truth Information - reco matched"<<endl;
		cout<<"Bare lepton"<<endl<<endl;
		for(Int_t i = 0; i <  (Int_t) lep_bareVec.size(); i++)
		{
			cout<<"Bare i: "<<i<<"\t";
			cout<<"Pt: "<< lep_bareVec[i].Pt()<<"\t";
			cout<<"Eta: "<< lep_bareVec[i].Eta()<<"\t";
			cout<<"Phi: "<< lep_bareVec[i].Phi()<<"\t";
			cout<<"M: "<< lep_bareVec[i].M()<<endl;
		}
		cout<<"Born lepton"<<endl;
		for(Int_t i = 0; i <  (Int_t) lep_bornVec.size(); i++)
		{
			cout<<"Born i: "<<i<<"\t";
			cout<<"Pt: "<< lep_bornVec[i].Pt()<<"\t";
			cout<<"Eta: "<< lep_bornVec[i].Eta()<<"\t";
			cout<<"Phi: "<< lep_bornVec[i].Phi()<<"\t";
			cout<<"M: "<< lep_bornVec[i].M()<<endl;
		}

		cout<<"mZ1_truth_bare: "<<higgs->mZ1_truth_bare<<endl;
		cout<<"mZ2_truth_bare: "<<higgs->mZ2_truth_bare<<endl;
		cout<<"m4l_truth_bare: "<<higgs->m4l_truth_bare<<endl;
		
		cout<<"mZ1_truth_born: "<<higgs->mZ1_truth_born<<endl;
		cout<<"mZ2_truth_born: "<<higgs->mZ2_truth_born<<endl;
		cout<<"m4l_truth_born: "<<higgs->m4l_truth_born<<endl;

		cout<<"mZ1_truth_dressed: "<<higgs->mZ1_truth_dressed<<endl;
		cout<<"mZ2_truth_dressed: "<<higgs->mZ2_truth_dressed<<endl;
		cout<<"m4l_truth_dressed: "<<higgs->m4l_truth_dressed<<endl;

		cout<<"Number of candidate FSR Photon found: "<<lepFSRTruthVec.size()<<endl;
		cout<<"# FSR combined with Z1: "<<ndRPassZ1<<endl;
		cout<<"# FSR combined with Z2: "<<ndRPassZ2<<endl;

		for(Int_t i = 0; i <  (Int_t) lepFSRTruthVec.size(); i++)
		{
			cout<<"FSR i: "<<i<<"\t";
			cout<<"Pt: "<< lepFSRTruthVec[i].Pt()<<"\t";
			cout<<"Eta: "<< lepFSRTruthVec[i].Eta()<<"\t";
			cout<<"Phi: "<< lepFSRTruthVec[i].Phi()<<endl;
		}

		cout<<"------------------"<<endl;

	}

	// Filling the trre turth info
	fillTruthInfo(higgs, bornSumVec);

	lep_bareVec.clear();
	lepBareTruthIndex.clear();
	lepBareTruthPDGID.clear();
	lepBornTruthIndex.clear();
	lepMotherTruthIndex.clear();
	lep_bornVec.clear();
	return;
}
// To Fill the true truth information
void HiggsAnalysis::fillTruthInfo(QuadLepton * higgs, TLorentzVector bornSumVec)
{
	// Truth PT
	Int_t higgsTruthIndex = -1;
	TLorentzVector higgsTruthVector;
	
	for (Int_t i = 0; i < (Int_t) event->mc.n(); i++)
	{
		Int_t pdgid = event->mc[i].pdgId();
		Int_t status = event->mc[i].status();
		if((pdgid==25 || pdgid == 39) && (status==2 || status==10902 || status==62)) 
		{
			higgsTruthIndex = i;
			higgsTruthVector.SetPtEtaPhiM(event->mc[i].pt(), 
							  			event->mc[i].eta(), 
							  			event->mc[i].phi(), 
							  			event->mc[i].m());
		}
	}

	// Fille the bornpt if no higgs has been found
	if(higgsTruthIndex != -1) higgs->truthVec = higgsTruthVector;
	else  higgs->truthVec = bornSumVec;
	
	vector<Int_t> trueBareIndex;
	vector<Int_t> trueBarePDGID;
	
	// Loop over the event to find born index
	for (Int_t i = 0; i < (Int_t) event->mc.n(); i++)
	{
		Int_t pdgid = event->mc[i].pdgId();
		Int_t status = event->mc[i].status();
		
		// Look for a Z that decays
		if(pdgid == 23 &&  ((status==22    || // Z->4l (Z coming from H->ZZ                     (l status 3 and 1)
                            status==155   || // Jimmy gg2ZZ NOT WORKING
                            status==10902 || // Pythia6 stable Z->ll (equivalent to status==2) (l status 1)
                            status==2     || // Pythia6 stable Z->ll                           (l status 1)
                            status==3     || // Pythia6 stable Z->ll                           (l status 3)
                            status==62)))    // PowHegPythia8 ZZ->4l                           (l status 3 and 1))
		{
			// To veto the third Z in ZH samples
			Bool_t higgsDecayZ = true;
			if(sampleProdType == sampleType::ZH) higgsDecayZ = checkParentHiggs(i);
			if(!higgsDecayZ) continue;

			for (Int_t j = 0; j < (Int_t) event->mc[i].child_index().size(); j++)
			{
				Int_t childIndex = event->mc[i].child_index().at(j);
				Int_t childpdgId = event->mc[childIndex].pdgId();
				Int_t childStatus = event->mc[childIndex].status();
				// If lepton, ask for bare right now
				if((abs(childpdgId) == 11 || abs(childpdgId) == 13) && childStatus == 1) 
				{
					trueBareIndex.push_back(childIndex);
					trueBarePDGID.push_back(childpdgId);
				}
				// If taus, ask for the one that decay's right now
				// Only impose Tau status 2 here...
				// This is like bare for Tau
				else if(abs(childpdgId) == 15 && (childStatus == 2 || childStatus == 10902))
				{
					trueBareIndex.push_back(childIndex);
					trueBarePDGID.push_back(childpdgId);
				}
			}
		}
	}

	// Sanity Check
	if(trueBareIndex.size() != 4) 
	{
		cout<<"Size trueBareIndex: "<<trueBareIndex.size()<<endl; 
		if(isDebugCall) printMCInfo(event->eventinfo.EventNumber());
		if(isDebugCall) for(Int_t i = 0; i <  (Int_t) trueBareIndex.size(); i++) cout<<i<<" trueBareIndex: "<<trueBareIndex[i]<<endl;
		return;
	}
	vector<Int_t> trueBornIndex = getBornIndexFromBare(trueBareIndex, trueBarePDGID, false);

	// Sanity Check
	if(trueBornIndex.size() != 4) {cout<<"Size trueBornIndex: "<<trueBornIndex.size()<<endl; return;}
	//for(Int_t i = 0; i < (Int_t) trueBornIndex.size(); i++) cout<<"index i: "<<i<<" Bare: "<<trueBareIndex[i]<<" Born: "<<trueBornIndex[i]<<endl;
	
	// Get the TLorentzVector
	vector<TLorentzVector> trueTurthBornVec;
	vector<Bool_t> usedBornVec; // to be used when reconstructing the Z
	
	for(Int_t i = 0; i < (Int_t) trueBornIndex.size(); i++) usedBornVec.push_back(false);


	vector<Int_t> trueComb1Index;
	vector<Int_t> trueComb2Index;

	trueComb1Index.push_back(trueBornIndex[0]);
	
	usedBornVec[0] = true;
	for(Int_t i = 0; i < (Int_t) trueBornIndex.size(); i++)
	{
		Int_t curr_Index = trueComb1Index[0];
		Int_t curr_pdgID =  event->mc[curr_Index].pdgId();
		// Don't want to reuse the vector
		if(usedBornVec[i]) continue;
		
		Int_t i_index = trueBornIndex[i];
		Int_t i_pdgID = event->mc[i_index].pdgId();

		// Check it is the same lepton and opposite charge
		// Easy way, sum of pdgid == 0
		if((i_pdgID + curr_pdgID ) != 0) continue;

		// Now ensure that the mother are the same
		// Check first if the size if equal
		if(event->mc[curr_Index].parent_index().size() != event->mc[i_index].parent_index().size()) continue;

		Bool_t sameParents = true;
		for(Int_t j = 0; j < (Int_t) event->mc[curr_Index].parent_index().size(); j++)
			if(event->mc[curr_Index].parent_index().at(j) !=  event->mc[i_index].parent_index().at(j)) sameParents = false;

		if(!sameParents) continue;

		// Now what is left is the right lepton
		trueComb1Index.push_back(i_index);
		
		usedBornVec[i] = true;
	}
	// The rest is comb2
	for(Int_t i = 0; i < (Int_t) trueBornIndex.size(); i++)
	{
		// Don't want to reuse the vector
		if(usedBornVec[i]) continue;
		trueComb2Index.push_back(trueBornIndex[i]);
		usedBornVec[i] = true;
	}

	// Sanity check
	if(trueComb1Index.size() != 2) cout<<"Comb1 size: "<<trueComb1Index.size() <<endl;
	if(trueComb2Index.size() != 2) cout<<"Comb2 size: "<<trueComb2Index.size() <<endl;

	// Reorder the leptons
	// First one is always lepPlus
	if(event->mc[trueComb1Index[0]].charge() == -1)
	{
		Int_t tempIndex = trueComb1Index[0];
		trueComb1Index[0] = trueComb1Index[1];
		trueComb1Index[1] = tempIndex;
	}
	if(event->mc[trueComb2Index[0]].charge() == -1)
	{
		Int_t tempIndex = trueComb2Index[0];
		trueComb2Index[0] = trueComb2Index[1];
		trueComb2Index[1] = tempIndex;
	}
	vector <TLorentzVector> trueComb1Vec = getVecFromIndex(trueComb1Index);
	vector <TLorentzVector> trueComb2Vec = getVecFromIndex(trueComb2Index);

	// Now find the closest to the pdgID mass
	Double_t diff1 = fabs((trueComb1Vec[0] + trueComb1Vec[1]).M() - pdgZMass);
	Double_t diff2 = fabs((trueComb2Vec[0] + trueComb2Vec[1]).M() - pdgZMass);
	
	if(trueComb1Vec.size() != 2) cout<<"trueComb1Vec size: "<<trueComb1Vec.size() <<endl;
	if(trueComb2Vec.size() != 2) cout<<"trueComb2Vec size: "<<trueComb2Vec.size() <<endl;

	vector<TLorentzVector> 	trueZ1Vec;
	vector<TLorentzVector> 	trueZ2Vec;

	if(diff1 < diff2)
	{
		trueZ1Vec = trueComb1Vec;
		trueZ2Vec = trueComb2Vec;
	}
	else
	{
		trueZ1Vec = trueComb2Vec;
		trueZ2Vec = trueComb1Vec;
	}


	higgs->getZ1()->getLepPlus()->m_momentumTruthTrueBorn	= trueZ1Vec[0];
	higgs->getZ1()->getLepNeg()->m_momentumTruthTrueBorn	= trueZ1Vec[1];
	higgs->getZ2()->getLepPlus()->m_momentumTruthTrueBorn	= trueZ2Vec[0];
	higgs->getZ2()->getLepNeg()->m_momentumTruthTrueBorn 	= trueZ2Vec[1];	

	if(isDebugCall)
	{
		cout<<"--------------------------------"<<endl;
		cout<<"Truth Information - true matched"<<endl;
		cout<<"Born lepton"<<endl<<endl;
		for(Int_t i = 0; i <  (Int_t) trueZ1Vec.size(); i++)
		{
			cout<<"Born i: "<<i<<"\t";
			cout<<"Pt: "<< trueZ1Vec[i].Pt()<<"\t";
			cout<<"Eta: "<< trueZ1Vec[i].Eta()<<"\t";
			cout<<"Phi: "<< trueZ1Vec[i].Phi()<<"\t";
			cout<<"M: "<< trueZ1Vec[i].M()<<endl;
		}
		for(Int_t i = 0; i <  (Int_t) trueZ2Vec.size(); i++)
		{
			cout<<"Born i: "<<i+2<<"\t";
			cout<<"Pt: "<< trueZ2Vec[i].Pt()<<"\t";
			cout<<"Eta: "<< trueZ2Vec[i].Eta()<<"\t";
			cout<<"Phi: "<< trueZ2Vec[i].Phi()<<"\t";
			cout<<"M: "<< trueZ2Vec[i].M()<<endl;
		}
		cout<<"--------------------------------"<<endl;
	}

	// Getting the PDGID to modify the branching ratio stuff
	vector<Int_t> trueBornPDGID;
	for(Int_t i = 0; i < (Int_t) trueBornIndex.size(); i++) 
		trueBornPDGID.push_back(abs(event->mc[trueBornIndex[i]].pdgId()));

	// now checking and modifying the branching ratio prodtype apporiately
	if((trueBornPDGID[0] == trueBornPDGID[1]) && (trueBornPDGID[1] == trueBornPDGID[2]) && (trueBornPDGID[2] == trueBornPDGID[3]) && (trueBornPDGID[3] == trueBornPDGID[0]))
	{
		if(higgs->getQuadType() == quadType::Mu4 || higgs->getQuadType() == quadType::Mu2El2) higgs->quadTypeBR = quadType::Mu4;
		if(higgs->getQuadType() == quadType::El4 || higgs->getQuadType() == quadType::El2Mu2) higgs->quadTypeBR = quadType::El4;	
	}
	else if((trueBornPDGID[0] !=  trueBornPDGID[3]) || (trueBornPDGID[1] != trueBornPDGID[2]))
	{
		if(higgs->getQuadType() == quadType::Mu4 || higgs->getQuadType() == quadType::Mu2El2) higgs->quadTypeBR = quadType::Mu2El2;
		if(higgs->getQuadType() == quadType::El4 || higgs->getQuadType() == quadType::El2Mu2) higgs->quadTypeBR = quadType::El2Mu2;	
	}
	else
		cout<<"Error: Branching ratio type cannot be decided from truth information"<<endl;
	// use the default if this can't be found....

}

// Checks if the parent is higgs or not...
Bool_t HiggsAnalysis::checkParentHiggs(Int_t index)
{
	if(event->mc[index].parent_index().size() == 0) return false; 

	for (Int_t i = 0; i < (Int_t )event->mc[index].parent_index().size(); i++)
	{
		Int_t currParentIndex = event->mc[index].parent_index().at(i);
		if( event->mc[currParentIndex].pdgId() == 25 ||  event->mc[currParentIndex].pdgId() == 39) return true;
	}

	return false;

}
// Recursively check if one of the parents in the chain is of the given PdgID
Bool_t HiggsAnalysis::checkParent(Int_t index, Int_t parentPDGID, Int_t constParentIndex)
{
	if(index == -1) return false;
	if(event->mc[index].parent_index().size() == 0) return false; 

	for (Int_t i = 0; i < (Int_t) event->mc[index].parent_index().size(); i++)
	{
		Int_t currParentIndex = event->mc[index].parent_index().at(i);
		if(currParentIndex == -1) continue;
		if(constParentIndex != -1 && constParentIndex != currParentIndex) continue;
		if(event->mc[currParentIndex].pdgId() ==  parentPDGID) return true;
	}

	Bool_t isUpTheChain = false;
	for (Int_t i = 0; i < (Int_t) event->mc[index].parent_index().size(); i++)
	{
		Int_t currParentIndex = event->mc[index].parent_index().at(i);
		isUpTheChain = isUpTheChain | checkParent(currParentIndex, parentPDGID, constParentIndex);
	}

	return isUpTheChain;
}

vector<TLorentzVector> HiggsAnalysis::getVecFromIndex(vector<Int_t> lepIndex)
{

	vector<TLorentzVector> lep_Vec;
	for(Int_t i = 0; i < (Int_t) lepIndex.size() && i < 4; i++)
	{
		Int_t currI = lepIndex[i];
		
		TLorentzVector l;
		if(currI != -1)
		{
			l.SetPtEtaPhiM(event->mc[currI].pt(), 
							  event->mc[currI].eta(), 
							  event->mc[currI].phi(), 
							  event->mc[currI].m());
		}
		lep_Vec.push_back(l);
	}
	return lep_Vec;
}

Int_t HiggsAnalysis::getIndexBarcodeMatch(Int_t truthBarcode, Int_t truthPDG)
{
	for(Int_t i = 0; i < event->mc.n(); i++)
	{
		if(event->mc[i].barcode() != truthBarcode) continue;

		// We have found the barcode lepton
		// Now make sure that it is a status 1 lepton with the same pdgID
		// If so these will be born leptons
		if(event->mc[i].status() == 1 && event->mc[i].pdgId() == truthPDG) return i;
		else cout<<"Mismatch in the MC container and truth level container"<<endl;

		break;
	}
	// For debugging purposes
	return -1;

}

Int_t HiggsAnalysis::getIndexDeltaRMatch(Int_t truthPDG, Double_t eta, Double_t phi)
{
	Double_t currDeltaR = 999;
	Double_t tempDeltaR = 999;
	Double_t currIndex = -1;

	for(Int_t i = 0; i < event->mc.n(); i++)
	{
		// Need status 1 leptons
		if(event->mc[i].status() != 1) continue;
		//cout<<"Here: "<<"TruthPDG: "<<truthPDG<<endl;
		// Need the same letpon
		if(fabs(event->mc[i].pdgId()) != fabs(truthPDG)) continue;

		// Finding the cloest one.
		tempDeltaR = DeltaR(eta, phi, event->mc[i].eta(), event->mc[i].phi());
		//cout<<"TempDeltaR: "<<tempDeltaR<<" index: "<<i<<endl;
		if(tempDeltaR < 0.05 )
		{
			if(tempDeltaR < currDeltaR)
			{
				currDeltaR = tempDeltaR;
				currIndex = i;
			}
		}
	}
	return currIndex;
}
vector<Int_t> HiggsAnalysis::getBornIndexFromBare(vector<Int_t> lepBareTruthIndex, vector<Int_t> lepBareTruthPDGID, Bool_t isRecoMatched)
{

	Bool_t isPythia6 = false;
	Bool_t isTauEvent = false;
	
	vector<Int_t> lepBornTruthIndex;
	// Copy the bare to born
	for(Int_t i = 0; i < (Int_t) lepBareTruthIndex.size(); i++) {lepBornTruthIndex.push_back(lepBareTruthIndex[i]);}

	// Now for every bare index, try to find a lepton with status 3, with has the same truth charge and the same mother
	for(Int_t i = 0; i < (Int_t) lepBornTruthIndex.size(); i++)
	{
		Bool_t isTauEventFilled = false;
		Int_t currI = lepBornTruthIndex[i];
		Int_t currPDGID = lepBareTruthPDGID[i];
		if(currI == -1) continue;

		Int_t currCharge = event->mc[currI].charge();
		
		if(event->mc[currI].parent_index().size() == 0) {continue;}
		if(event->mc[currI].parent_index().size() != 1) 
			{cout<<"Error fillTruthRecoMatchedInfo: Parent Size is "<<event->mc[currI].parent_index().size() <<" . Choosing frist parent. Curr Index: "<<currI<<endl;}

		Int_t currParentIndex = event->mc[currI].parent_index().at(0);
		// If parent is a tau, loop over to find a same charge tau with status 3
		if(abs(event->mc[currParentIndex].pdgId()) == 15)
		{
			lepBornTruthIndex[i] = currParentIndex;
			currI = currParentIndex;
			currPDGID = event->mc[currI].pdgId();
			currCharge = event->mc[currI].charge();
			if(event->mc[currI].parent_index().size() == 0) {continue;}
			if(event->mc[currI].parent_index().size() != 1) 
				{cout<<"Error fillTruthRecoMatchedInfo: Parent Size is "<<event->mc[currI].parent_index().size() <<" . Choosing frist parent. Curr Index: "<<currI<<endl;}

			currParentIndex = event->mc[currI].parent_index().at(0);
			isTauEvent = true;
		}
		
		if(isTauEventFilled) continue;	

		// Now looping over the mc Container to find the born lepton
		for(Int_t j = 0; j < event->mc.n(); j++)
		{
			// First compare just the charge, status and PDGID
			if(event->mc[j].pdgId() == currPDGID &&
			   event->mc[j].charge() == currCharge &&
			   event->mc[j].status() == 3)
			{
				// For pythia6, born and bare don't share the same mother
				// So look for a mother that is a Z with status 3
				for(Int_t k = 0; k < (Int_t) event->mc[j].parent_index().size(); k++)
				{
					Int_t motherIndex = event->mc[j].parent_index().at(k);
					if(event->mc[motherIndex].pdgId() == 23 && event->mc[motherIndex].status() == 3) 
					{
						// Check if the index has been used already
						if(std::find(lepBornTruthIndex.begin(), lepBornTruthIndex.end(),j) != lepBornTruthIndex.end()) continue;
						isPythia6 = true;
					}
				}
				// This works for pythia8
				// Now check the mother
				for(Int_t k = 0; k < (Int_t) event->mc[j].parent_index().size(); k++)
				{	
					if(event->mc[j].parent_index().at(k) == currParentIndex) 
					{
						lepBornTruthIndex[i] = j; 
						break;
					}
				}
			}
		}
	}

	// Temporary fix for pythia 6
	if(isPythia6)
	{
		lepBornTruthIndex = getBornIndexFromBarePythia6(lepBareTruthIndex, lepBareTruthPDGID, isRecoMatched);
	}
	return lepBornTruthIndex;
}
// Special functions for pythia 6
vector<Int_t> HiggsAnalysis::getBornIndexFromBarePythia6(vector<Int_t> lepBareTruthIndex, vector<Int_t> lepBareTruthPDGID, Bool_t isRecoMatched)
{
	vector<Int_t> trueBareIndex;
	vector<Int_t> trueBarePDGID;
	vector<Int_t> trueBornIndex;
	vector<Bool_t> isTauMother;
	Bool_t isTauEvent = false;
	// Loop over the event to find born index
	for (Int_t i = 0; i < (Int_t) event->mc.n(); i++)
	{
		Int_t pdgid = event->mc[i].pdgId();
		Int_t status = event->mc[i].status();
		
		// Look for a Z that decays
		if(pdgid == 23 &&  ((status==22    || // Z->4l (Z coming from H->ZZ                     (l status 3 and 1)
                            status==155   || // Jimmy gg2ZZ NOT WORKING
                            status==10902 || // Pythia6 stable Z->ll (equivalent to status==2) (l status 1)
                            status==2     || // Pythia6 stable Z->ll                           (l status 1)
                            status==3     || // Pythia6 stable Z->ll                           (l status 3)
                            status==62)))    // PowHegPythia8 ZZ->4l                           (l status 3 and 1))
		{
			for (Int_t j = 0; j < (Int_t) event->mc[i].child_index().size(); j++)
			{
				Int_t childIndex = event->mc[i].child_index().at(j);
				Int_t childpdgId = event->mc[childIndex].pdgId();
				Int_t childStatus = event->mc[childIndex].status();
				// If lepton, ask for bare right now
				if((abs(childpdgId) == 11 || abs(childpdgId) == 13) && childStatus == 1) 
				{
					trueBareIndex.push_back(childIndex);
					trueBarePDGID.push_back(childpdgId);
					isTauMother.push_back(false);

					// Copy the index into born for now as well
					trueBornIndex.push_back(childIndex);
					
				}
				// If taus, ask for the one that decay's right now
				else if(abs(childpdgId) == 15 && childStatus == 2) //status efore decay
				{
					if(!isRecoMatched)
					{
						isTauEvent = true;
						trueBareIndex.push_back(childIndex);
						trueBarePDGID.push_back(childpdgId);
						isTauMother.push_back(false);

						// Copy the index into born for now as well
						trueBornIndex.push_back(childIndex);
					}
					else
					{
						// Protection against a stupid sample where the tau didn't decay right....
						Bool_t foundLeptoninTauDecay = false;
						for (Int_t k = 0; k < (Int_t) event->mc[childIndex].child_index().size(); k++)
						{
							Int_t childIndexTau = event->mc[childIndex].child_index().at(k);
							Int_t childpdgIdTau = event->mc[childIndexTau].pdgId();
							Int_t childStatusTau = event->mc[childIndexTau].status();
							// If lepton, ask for bare right now
							if((abs(childpdgIdTau) == 11 || abs(childpdgIdTau) == 13) && childStatusTau == 1 ) 
							{
								foundLeptoninTauDecay = true;
								trueBareIndex.push_back(childIndexTau);
								trueBarePDGID.push_back(childpdgIdTau);
								isTauMother.push_back(true);
								// Copy the index into born for now as well
								trueBornIndex.push_back(childIndexTau);
							}

						}
						if(!foundLeptoninTauDecay)
						{
							trueBareIndex.push_back(-1);
							trueBarePDGID.push_back(-1);
							isTauMother.push_back(false);
							// Copy the index into born for now as well
							trueBornIndex.push_back(-1);

						}
					}
				}
			}
		}
	}

	// Now loopover the born
	// if not a tau mother, find the sequential born lepton
	for(Int_t k = 0; k <  (Int_t) trueBornIndex.size(); k++)
	{
		Bool_t filledTruth = false;	
		if(trueBornIndex[k] == -1) continue;
		// If tau mother replace the born index with taus (status 3, if it exits)
		if(isTauMother[k]) 
		{
			Int_t currI = trueBornIndex[k];
			for (Int_t i = 0; i < (Int_t) event->mc.n(); i++)
			{
				Int_t pdgid = abs(event->mc[i].pdgId());
				Int_t status = event->mc[i].status();
				if(pdgid != 15) continue;

				//if((hasStatus3Taus && status == 3) || (!hasStatus3Taus && status == 2))
				if(status == 2)
				{
					if(event->mc[i].charge() == event->mc[currI].charge())
					{ 
						trueBornIndex[k] = i;
						break;
					}
				}
			}
		 //continue;
		}
		Int_t currPDGID = event->mc[trueBornIndex[k]].pdgId();
		// Loop over the event to find born index
		for (Int_t i = 0; i < (Int_t) event->mc.n(); i++)
		{
			Int_t pdgid = event->mc[i].pdgId();
			Int_t status = event->mc[i].status();
			
			// Look for a Z that decays
			if(pdgid == 23 &&  ((status==22    || // Z->4l (Z coming from H->ZZ                     (l status 3 and 1)
                            status==155   || // Jimmy gg2ZZ NOT WORKING
                            status==10902 || // Pythia6 stable Z->ll (equivalent to status==2) (l status 1)
                            status==2     || // Pythia6 stable Z->ll                           (l status 1)
                            status==3     || // Pythia6 stable Z->ll                           (l status 3)
                            status==62)))    // PowHegPythia8 ZZ->4l                           (l status 3 and 1))
			{
				for (Int_t j = 0; j < (Int_t) event->mc[i].child_index().size(); j++)
				{
					Int_t childIndex = event->mc[i].child_index().at(j);
					Int_t childpdgId = event->mc[childIndex].pdgId();
					Int_t childStatus = event->mc[childIndex].status();

					if(abs(currPDGID) != 15 && (currPDGID != childpdgId)) continue;
					// If lepton, ask for born right now
					//if((abs(childpdgId) == 11 || abs(childpdgId) == 13 || abs(childpdgId) == 15) && childStatus == 3) 
					// 
					if((abs(childpdgId) == 11 || abs(childpdgId) == 13 || abs(childpdgId) == 15) && childStatus == 3) 
					{
						if(std::find(trueBornIndex.begin(), trueBornIndex.end(), childIndex) != trueBornIndex.end()) continue;
						trueBornIndex[k] = childIndex;
						filledTruth = true;
						break;
					}
				}
			}
			if(filledTruth)break;
		}
	}

	if(trueBornIndex.size() < 4)
	{
		for(Int_t i = (Int_t) trueBornIndex.size(); i < 4; i++)
		{
			trueBornIndex.push_back(-1);
			trueBareIndex.push_back(-1);
			
		}
	}
	vector<Int_t> lepBornTruthIndex;
	Int_t lepBornTruthIndexTemp [4];
	for(Int_t i = 0; i < (Int_t) lepBareTruthIndex.size() && i < 4; i++)
	{
		// Default value
		lepBornTruthIndexTemp[i] = lepBareTruthIndex[i];
		for(Int_t j = 0; j <  (Int_t) trueBareIndex.size(); j++)
		{
			if(trueBareIndex[j] == lepBareTruthIndex[i])
			{
				// push back the corresponding born index
				lepBornTruthIndexTemp[i] = trueBornIndex[j];
				break;
			}	
		}
	}
	if(trueBareIndex.size() != 4) 
	{
		//lepBornTruthIndex = lepBareTruthIndex;
		//return lepBornTruthIndex;
	}
	for(Int_t i = 0; i < 4; i++)
	{
		lepBornTruthIndex.push_back(lepBornTruthIndexTemp[i]);
	//	cout<<"lepBornTruthIndex "<<i<<" :"<< lepBornTruthIndexTemp[i]<<endl;
	}
	if(isDebugCall)
	{
		//printMCInfo(event->eventinfo.EventNumber());
		for(Int_t i = 0; i < (Int_t) trueBareIndex.size(); i++)
		{
			cout<<i<<" trueBareIndex :"<< trueBareIndex[i]<<" trueBornIndex: ";
			cout<< trueBornIndex[i]<<" lepBareTruthIndex: "<<lepBareTruthIndex[i];
			cout<< " lepBornTruthIndex: "<<lepBornTruthIndex[i]<<endl;
		}
		cout<<endl;
	}


	if(lepBornTruthIndex.size() != 4 || lepBareTruthIndex.size() != 4) 
	{
		cout<<"Pythia 6: Event Number: "<<event->eventinfo.EventNumber()<<" Size of the lepBornTruthIndex: "<<lepBornTruthIndex.size() <<endl;
		cout<<"Pythia 6: Event Number: "<<event->eventinfo.EventNumber()<<" Size of the lepBareTruthIndex: "<<lepBareTruthIndex.size() <<endl;
		
	}
	return lepBornTruthIndex;

}

vector<Int_t> HiggsAnalysis::getMotherIndexFromBare(vector<Int_t> lepBareTruthIndex, vector<Int_t> lepBareTruthPDGID, vector<Int_t> lepBornTruthIndex)
{
	//vector<Int_t> lepBornTruthIndex;
	vector<Int_t> lepMotherTruthIndex;

	// Copy the bare to born
	//for(Int_t i = 0; i < (Int_t) lepBareTruthIndex.size(); i++) lepBornTruthIndex.push_back(lepBareTruthIndex[i]);

	// Now for every bare index, try to find a lepton with status 3, with has the same truth charge and the same mother
	for(Int_t i = 0; i < (Int_t) lepBornTruthIndex.size(); i++)
	{
		Int_t currI = lepBornTruthIndex[i];
		if(currI == -1) continue;
		
		if(abs(event->mc[currI].pdgId()) == 15){lepMotherTruthIndex.push_back(currI);}
		if(event->mc[currI].parent_index().size() == 0) {continue;}
		if(event->mc[currI].parent_index().size() != 1) 
			{cout<<"Error fillTruthRecoMatchedInfo: Parent Size is "<<event->mc[currI].parent_index().size() <<" . Choosing frist parent. Curr Index: "<<currI<<endl;}

		Int_t currParentIndex = event->mc[currI].parent_index().at(0);
		if(abs(event->mc[currParentIndex].pdgId()) == 15) 
		{
			if(event->mc[currParentIndex].parent_index().size() != 0) 
			{		
				Int_t currParentIndexTau = event->mc[currParentIndex].parent_index().at(0);
				lepMotherTruthIndex.push_back(currParentIndexTau);	
			}
		}
		lepMotherTruthIndex.push_back(currParentIndex);
	}
	for(Int_t i = 0; i < (Int_t) lepBareTruthIndex.size(); i++)
	{
		Int_t currI = lepBareTruthIndex[i];
		if(currI == -1) continue;
		
		if(abs(event->mc[currI].pdgId()) == 15){lepMotherTruthIndex.push_back(currI);}
		if(event->mc[currI].parent_index().size() == 0) {continue;}
		if(event->mc[currI].parent_index().size() != 1) 
			{cout<<"Error fillTruthRecoMatchedInfo: Parent Size is "<<event->mc[currI].parent_index().size() <<" . Choosing frist parent. Curr Index: "<<currI<<endl;}

		Int_t currParentIndex = event->mc[currI].parent_index().at(0);
		if(abs(event->mc[currParentIndex].pdgId()) == 15) 
		{	
			if(event->mc[currParentIndex].parent_index().size() != 0) 
			{		
				Int_t currParentIndexTau = event->mc[currParentIndex].parent_index().at(0);
				lepMotherTruthIndex.push_back(currParentIndexTau);	
			}
		}
		lepMotherTruthIndex.push_back(currParentIndex);
	}

	//if(isDebugCall)
	//{
	//	for(Int_t i = 0; i < (Int_t) lepMotherTruthIndex.size(); i++) cout<<"lepMotherTruthIndex: "<<lepMotherTruthIndex[i]<<endl;
	//}


	return lepMotherTruthIndex;
}
// Assumption index1 is for the fsr Photon and index 2 is for the lepton
Bool_t HiggsAnalysis::checkSameParent(Int_t index1, Int_t index2)
{
	Bool_t sameParents = true;
	if(index1 == -1 || index2 == -1) return false;

	if(abs(event->mc[index1].pdgId()) == 15) index1 = event->mc[index1].parent_index().at(0);
	if(abs(event->mc[index2].pdgId()) == 15) index2 = event->mc[index2].parent_index().at(0);

	vector <Int_t> leptonParent;

	// Push the parent in the vector
	for(Int_t j = 0; j < (Int_t) event->mc[index2].parent_index().size(); j++)
	{
		Int_t parentIndex = event->mc[index2].parent_index().at(j);
		leptonParent.push_back(parentIndex);
		if(abs(event->mc[parentIndex].pdgId()) == 15)
		{
			for(Int_t k = 0; k < (Int_t) event->mc[parentIndex].parent_index().size(); k++)
			{
				Int_t parentParentIndex = event->mc[parentIndex].parent_index().at(k);
				leptonParent.push_back(parentParentIndex);
			}

		}
			
	}
	//if(isDebugCall)
	//{
	//	cout<<"Parents: \t";
	//	for (Int_t i = 0; i < (Int_t) leptonParent.size(); i++) cout<<leptonParent[i]<<"\t";
	//	cout<<endl;

	//	cout<<"Photon Parents: \t";
	//	for (Int_t i = 0; i <  (Int_t) event->mc[index1].parent_index().size(); i++) cout<<event->mc[index1].parent_index().at(i)<<"\t";
	//	cout<<endl;

	//}

	for(Int_t j = 0; j < (Int_t) event->mc[index1].parent_index().size(); j++)
	{
		if(std::find(leptonParent.begin(), leptonParent.end(), event->mc[index1].parent_index().at(j)) != leptonParent.end()) continue;
		return false;
	}

	return true;
}

void HiggsAnalysis::fillTruthJetsInfo(QuadLepton * higgs)
{
	ChargedLepton *leadingJet = 0;
	Bool_t filledLeadingJet = false;
	// For Jet Information

	Double_t leadingJetPt = -999*1000;
	for(Int_t i = 0; i < (Int_t) jetsTruthEvent.size(); i++)
	{
		// Check leading jet
		if(jetsTruthEvent[i]->get4Momentum()->Pt() > leadingJetPt)
		{
			// save the leading jet
			leadingJet = jetsTruthEvent[i];
			leadingJetPt = jetsTruthEvent[i]->get4Momentum()->Pt();
			filledLeadingJet = true;
		}
	}

	higgs->n_jets_truth_bare = jetsTruthEvent.size();
	higgs->leading_jet_pt_truth_bare = leadingJetPt/1000;

	if(isDebugCall)
	{
		cout<<"------------------------"<<endl;
		cout<<"Truth Jet: "<<higgs->n_jets_truth_bare<<endl;		
		cout<<"Truth Leading Jet pT: "<<higgs->leading_jet_pt_truth_bare<<endl;
		cout<<"------------------------"<<endl;
		
	}

}
// Get information on the sample
// see if it has a status 3 taus
Bool_t HiggsAnalysis::containStatus3Taus()
{
	for (Int_t i = 0; i < (Int_t) event->mc.n(); i++)
	{
		Int_t pdgid = event->mc[i].pdgId();
		Int_t status = event->mc[i].status();
		
		// If taus, ask for the one that decay's right now
		if(abs(pdgid) == 15 && status == 3)
		{
			return true;
		}
	}
	return false;

}

// Get the truth level quad type:
Int_t HiggsAnalysis::getTruthQuadType()
{
	Int_t nTau = 0;
	Int_t nMu = 0;
	Int_t nEl = 0;

	if(!isMC) return truthTypeQuad::noStatus;
	// Loop over the event to find born index
	for (Int_t i = 0; i < (Int_t) event->mc.n(); i++)
	{
		Int_t pdgid = event->mc[i].pdgId();
		Int_t status = event->mc[i].status();
		
		// Look for a Z that decays
		if(pdgid == 23 &&  ((status==22    || // Z->4l (Z coming from H->ZZ                     (l status 3 and 1)
                            status==155   || // Jimmy gg2ZZ NOT WORKING
                            status==10902 || // Pythia6 stable Z->ll (equivalent to status==2) (l status 1)
                            status==2     || // Pythia6 stable Z->ll                           (l status 1)
                            status==3     || // Pythia6 stable Z->ll                           (l status 3)
                            status==62)))    // PowHegPythia8 ZZ->4l                           (l status 3 and 1))
		{
			// To veto the third Z in ZH samples
			Bool_t higgsDecayZ = true;
			if(sampleProdType == sampleType::ZH) higgsDecayZ = checkParentHiggs(i);
			if(!higgsDecayZ) continue;

			for (Int_t j = 0; j < (Int_t) event->mc[i].child_index().size(); j++)
			{
				Int_t childIndex = event->mc[i].child_index().at(j);
				Int_t childpdgId = event->mc[childIndex].pdgId();
				Int_t childStatus = event->mc[childIndex].status();

				// check is parent is higgs
				if(!checkParent(childIndex, 25)) continue;
				
				// If lepton, ask for bare right now
				if(abs(childpdgId) == 11  && childStatus == 1) 
				{
					nEl ++;
				}
				if(abs(childpdgId) == 13 && childStatus == 1) 
				{
					nMu++;
				}
				// If taus, ask for the one that decay's right now
				// Only impose Tau status 2 here...
				// This is like bare for Tau
				else if(abs(childpdgId) == 15 && (childStatus == 2 || childStatus == 10902))
				{
					nTau++;
				}
			}
		}
	}

	if(nTau >= 4) return truthTypeQuad::_4tau;
	else if(nEl >= 4) return truthTypeQuad::_4e;
	else if(nMu >= 4) return truthTypeQuad::_4mu;
	else if (nTau >= 2 && nMu >= 2) return truthTypeQuad::_2tau2mu;
	else if (nTau >= 2 && nEl >= 2) return truthTypeQuad::_2tau2e;
	else if (nMu >= 2 && nEl >= 2) return truthTypeQuad::_2mu2e;
	else return truthTypeQuad::noStatus;

}

// Helper function to fill lepID
void HiggsAnalysis::fillLepID(QuadLepton * higgs)
{
	fillLepIDLepton(higgs->getZ1()->getLepPlus());
	fillLepIDLepton(higgs->getZ1()->getLepNeg());
	fillLepIDLepton(higgs->getZ2()->getLepPlus());
	fillLepIDLepton(higgs->getZ2()->getLepNeg());
}
// Filles in LepID
void HiggsAnalysis::fillLepIDLepton(ChargedLepton * lep)
{
	if(lep->getFlavor() == flavor::Electron)
	{
		if(dataYear == 2011)
		{
			lep->lepID = leptonIDType::el_loosepp_H4l_Ep_comb;

			// if not EP change the type to the default on
			if(fabs(lep->GetElectron()->cl_pt() - lep->clPt) < 0.00001) 
				lep->lepID = leptonIDType::el_loosepp_H4l;
		}
		else if(dataYear == 2012)
		{
			if(useLikelihood) 
			{
				lep->lepID = leptonIDType::el_likelihood_loose_Ep_comb;

				// if not EP change the type to the default on
				if(fabs(lep->GetElectron()->cl_pt() - lep->clPt) < 0.00001) 
					lep->lepID = leptonIDType::el_likelihood_loose;

			}
			else 
			{
				lep->lepID = leptonIDType::el_multilepton_Ep_comb;
				// if not EP change the type to the default on
				if(fabs(lep->GetElectron()->cl_pt() - lep->clPt) < 0.00001) 
					lep->lepID = leptonIDType::el_multilepton;

			}
		}
	}
	else if(lep->getFlavor() == flavor::Muon)
	{
		if(lep->getType() == leptonType::MuonCalo) lep->lepID = leptonIDType::mu_calomuon;
		else
		{
			D3PDReader::MuonD3PDObjectElement* mu =  lep->GetMuon();

			if(mu->isCombinedMuon()) lep->lepID = leptonIDType::mu_staco_cb;
			else if(mu->isSegmentTaggedMuon())	lep->lepID = leptonIDType::mu_staco_st;			
			else if(mu->isStandAloneMuon()) lep->lepID = leptonIDType::mu_staco_sa;
		}
	}
}

// Muon ID and ms vars
void HiggsAnalysis::fillMuonHelperVars()
{
	D3PDReader::MuonD3PDObject * mu = &(event->mu_staco);
	for(Int_t i = 0; i < mu->n(); i++)
	{
		Int_t isStandAloneMu = (*mu)[i].isStandAloneMuon();
		
		Double_t ptMs = (1/fabs((*mu)[i].me_qoverp())*sin((*mu)[i].me_theta()));
		Double_t phiMs = (*mu)[i].me_phi();
		Double_t etaMs = -TMath::Log(tan((*mu)[i].me_theta()/2));
		Double_t pTId;
		Double_t phiId;
		Double_t etaId;
		if(isStandAloneMu) 
		{
			pTId = (*mu)[i].pt();
			phiId = (*mu)[i].phi();
			etaId = (*mu)[i].eta();
		}
		else{
			pTId = (1/fabs((*mu)[i].id_qoverp())*sin((*mu)[i].id_theta())) ;
			phiId = (*mu)[i].id_phi();
			etaId = -TMath::Log(tan((*mu)[i].id_theta()/2));
		}
		
		// For segment tagged
		if((*mu)[i].isSegmentTaggedMuon()) 
		{
			ptMs = (*mu)[i].pt();
			phiMs = (*mu)[i].phi();
			etaMs = (*mu)[i].eta();
		}

		(*mu)[i].me_pt  = ptMs;
		(*mu)[i].me_phi() = phiMs;
		(*mu)[i].me_eta = etaMs;

		(*mu)[i].id_pt  = pTId;
		(*mu)[i].id_phi() = phiId;
		(*mu)[i].id_eta = etaId;

		(*mu)[i].id_pt_unsmeared = (*mu)[i].id_pt;
		(*mu)[i].me_pt_unsmeared = (*mu)[i].me_pt;
		(*mu)[i].cb_pt_unsmeared = (*mu)[i].pt();
		//cout<<"Event "<<event->eventinfo.EventNumber()<<endl;
		//cout<<"staco bfID theta: "<<(*mu)[i].id_theta()/2<<" MS theta: "<<(*mu)[i].me_theta()/2<<endl;
		//cout<<"staco bfID eta: "<<(*mu)[i].id_eta<<" MS eta: "<<(*mu)[i].me_eta<<endl;
		//
		if(isDebugCall)
		{
			cout<<"Staco i: CB pt: "<<(*mu)[i].cb_pt_unsmeared<<" ID: "<<(*mu)[i].id_pt_unsmeared<<" MS: "<<(*mu)[i].me_pt_unsmeared<<endl;
		}
	}

	mu = &(event->mu_calo);
	for(Int_t i = 0; i < mu->n(); i++)
	{
		
		Double_t ptMs = (*mu)[i].pt();
		Double_t phiMs = (*mu)[i].phi();
		Double_t etaMs = (*mu)[i].eta();
		Double_t pTId;
		Double_t phiId;
		Double_t etaId;
		
		pTId = (1/fabs((*mu)[i].id_qoverp())*sin((*mu)[i].id_theta())) ;
		phiId = (*mu)[i].id_phi();
		etaId = -TMath::Log(tan((*mu)[i].id_theta()/2));

		(*mu)[i].me_pt  = ptMs;
		(*mu)[i].me_phi() = phiMs;
		(*mu)[i].me_eta = etaMs;

		(*mu)[i].id_pt  = pTId;
		(*mu)[i].id_phi() = phiId;
		(*mu)[i].id_eta = etaId;

		(*mu)[i].id_pt_unsmeared = (*mu)[i].id_pt;
		(*mu)[i].me_pt_unsmeared = (*mu)[i].me_pt;
		(*mu)[i].cb_pt_unsmeared = (*mu)[i].pt();

		//cout<<"calo  bfID theta: "<<(*mu)[i].id_theta()/2<<" MS theta: "<<(*mu)[i].me_theta()/2<<endl;
		//cout<<"calo bfID eta: "<<(*mu)[i].id_eta<<" MS eta: "<<(*mu)[i].me_eta<<endl;
		if(isDebugCall)
		{
			cout<<"Calo i: CB pt: "<<(*mu)[i].cb_pt_unsmeared<<" ID: "<<(*mu)[i].id_pt_unsmeared<<" MS: "<<(*mu)[i].me_pt_unsmeared<<endl;
		}

		
	}

	for(Int_t i = 0; i < el_cur->n(); i++)
	{
		
		Double_t eta_trk = (&((*el_cur)[i]))->tracketa();
		Double_t E = (&((*el_cur)[i]))->cl_E();
		Double_t pT = E/cosh(eta_trk);	

		// Just storage
		(&((*el_cur)[i]))->pT_unsmeared = pT;
	}
}

////////////////////////////////////////////////////////////////////////////////////////
//				To Clean the vars that are event Dependant...
////////////////////////////////////////////////////////////////////////////////////////
void HiggsAnalysis::clearVar()
{
	// Clean the Vector that contains the muons
	//while(!muOverlap.empty()) delete muOverlap.back(), muOverlap.pop_back();
	//while(!elOverlap.empty()) delete elOverlap.back(), elOverlap.pop_back();
	//while(!elLooseOverlap.empty()) delete elLooseOverlap.back(), elLooseOverlap.pop_back();
	//while(!jetsOverlap.empty()) delete jetsOverlap.back(), jetsOverlap.pop_back();
	//while(!jetsTruthEvent.empty()) delete jetsTruthEvent.back(), jetsTruthEvent.pop_back();
	//while(!jetsOverlap_Fid.empty()) delete jetsOverlap_Fid.back(), jetsOverlap_Fid.pop_back();
	//while(!jetsTruthEvent_Fid.empty()) delete jetsTruthEvent_Fid.back(), jetsTruthEvent_Fid.pop_back();
	while(!eventToClean4mu.empty()) delete eventToClean4mu.back(), eventToClean4mu.pop_back();
	while(!eventToClean4el.empty()) delete eventToClean4el.back(), eventToClean4el.pop_back();
	while(!eventToClean2mu2el.empty()) delete eventToClean2mu2el.back(), eventToClean2mu2el.pop_back();
	while(!eventToCleanCR.empty()) delete eventToCleanCR.back(), eventToCleanCR.pop_back();
	while(!eventToCleanDi.empty()) delete eventToCleanDi.back(), eventToCleanDi.pop_back();
	muEvent.clear();
	elEvent.clear();
	elLooseEvent.clear();
	jetsEvent.clear();
	muOverlap.clear();
	elOverlap.clear();
	elLooseOverlap.clear();
	jetsOverlap.clear();
	jetsTruthEvent.clear();
	jetsEvent_Fid.clear();
	jetsTruthEvent_Fid.clear();
	jetsOverlap_Fid.clear();
	eventToClean4mu.clear();
	eventToClean4el.clear();
	eventToClean2mu2el.clear();
	eventToCleanCR.clear();
}

// Helper var
Double_t HiggsAnalysis::DeltaR (Double_t eta_1, Double_t phi_1, Double_t eta_2, Double_t phi_2)
{
	Double_t dR=0;
	Double_t eta2 = (eta_1-eta_2)*(eta_1-eta_2);
	Double_t tmp_dphi = (fabs(phi_1-phi_2) > TMath::Pi()) ? 2*TMath::Pi()-fabs(phi_1-phi_2) : fabs(phi_1-phi_2);
	Double_t phi2 = tmp_dphi*tmp_dphi;
	dR = sqrt( eta2 + phi2 );
	return dR;
}

// Gets the Delta R for calo Iso purposes
Double_t HiggsAnalysis::DeltaRCaloIso (ChargedLepton *lep1, ChargedLepton *lep2)
{
	// Just to Make life a bit easier
	ChargedLepton *curr_lep [2];
	curr_lep[0] = lep1;
	curr_lep[1] = lep2;
	// Varible to compute the DeltaR
	Double_t eta[2];
	Double_t phi[2];

	for(Int_t i = 0; i< 2; i++)
	{
		if(curr_lep[i]->getFlavor() == flavor::Electron)
		{
			//if(dataYear == 2011)
			//{
			//	eta[i] = curr_lep[i]->GetElectron()->tracketa();
			//	phi[i] = curr_lep[i]->GetElectron()->trackphi();	
			//}
			//else if (dataYear == 2012)
			//{
			//	eta[i] = curr_lep[i]->GetElectron()->etas2();
			//	phi[i] = curr_lep[i]->GetElectron()->phis2();
			//}
			
			eta[i] = curr_lep[i]->GetElectron()->etas2();
			phi[i] = curr_lep[i]->GetElectron()->phis2();
		}
		else if(curr_lep[i]->getFlavor() == flavor::Muon)
		{
				eta[i] = curr_lep[i]->get4Momentum()->Eta();
				phi[i] = curr_lep[i]->get4Momentum()->Phi();
		}
	}

	return DeltaR(eta[0], phi[0], eta[1], phi[1]);
}
// This function returns the MC higgs mass that is in the file name
// If it is background, it returns -1
Double_t HiggsAnalysis::getMCHiggsMass()
{
	// Storing the name of the file
	TString name;
	name = currFileName;
	// Overwrite from grid
	if(runningGrid) name = gridFileName;

	// Check if no tau sample
	if(name.Contains("noTau")) noTauSample = true;
	else noTauSample = false;

	// Check what type of Sample it is
	// Check if no tau sample
	if(name.Contains("Pythia") ) generatorName = MCGeneratorName::Pythia;

		// This is for JHU samples
	Int_t RunNumber = 0;
	if(isMC) RunNumber = event->eventinfo.mc_channel_number();
	if((RunNumber>=169716 && RunNumber<=169717) ||  // 2011 JHU
       (RunNumber>=167604 && RunNumber<=167605) ||  // 2011 JHU
       (RunNumber>=167607 && RunNumber<=167607) ||  // 2012 JHU
       (RunNumber>=169710 && RunNumber<=169711) ||  // 2012 JHU
       (RunNumber>=167124 && RunNumber<=167125) ||  // 2012 JHU
       (RunNumber>=167127 && RunNumber<=167127) ||
	   (RunNumber>=167600 && RunNumber<=167603) || // 2011 JHU
       (RunNumber>=167606 && RunNumber<=167606) || // 2011 JHU
       (RunNumber>=167120 && RunNumber<=167123) || // 2012 JHU
       (RunNumber>=167126 && RunNumber<=167126) ) 
	{  		  
		noTauSample = true;
		generatorName = MCGeneratorName::Pythia;
	} 

	// Splitting the file path
	TObjArray *parts = name.Tokenize(".");
	vector<TString> partName;
	if(parts->GetEntriesFast()) {
	   TIter iString(parts);
	   TObjString* os=0;
	   while ((os=(TObjString*)iString())) {
	  	partName.push_back(os->GetString());
	   }
	}
	
	// To loop over all the parts
	Int_t nMax = partName.size();

	for(Int_t i = 0; i < nMax; i++)
	{
		// Getting the first part of the container name
		if(partName[i].Contains("mc11_7TeV") || partName[i].Contains("mc12_8TeV"))
		{
			//cout << name << ":  ";
			// Getting the mc Run Number
			Int_t j = i + 1;
			TString numSample = partName[j];
			mcRunNumber = numSample.Atoi();

			// Higgs mass is 2 aways
			Int_t indexStart = 0;
			i = i+2;
			// Find the sample type and the start index of the mass
			
			// To take care of the madgraph z'z' samples
			if(partName[i].Contains("_ZZp") || partName[i].Contains("_ZpZp") )
			{
				indexStart = partName[i].Index("ggH") + 3;				
				sampleProdType = sampleType::ggF_ZpZp;
				if(curEvent == 0)
				{
					cout<<"---------------------------"<<endl;
					cout<<"This is a Z(')Z' sample"<<endl;
					cout<<"---------------------------"<<endl;
				}

			}
			
			else if(partName[i].Contains("ggH"))
			{
				//cout<<"ggH Sample ";
				indexStart = partName[i].Index("ggH") + 3;
				sampleProdType = sampleType::ggF;
				if(curEvent == 0)
				{
					cout<<"---------------------------"<<endl;
					cout<<"This is a ggF sample"<<endl;
					cout<<"---------------------------"<<endl;
				}
			}
			else if(partName[i].Contains("VBF"))
			{
				indexStart = partName[i].Index("VBF") + 3;

				// workaround for VBFH sample
				if(partName[i].Contains("VBFH")) indexStart++;
				sampleProdType = sampleType::VBF;
				if(curEvent == 0)
				{
					cout<<"---------------------------"<<endl;
					cout<<"This is a VBF sample"<<endl;
					cout<<"---------------------------"<<endl;
				}
			}

			else if(partName[i].Contains("WH"))
			{
				//cout<<"WH Sample ";
				indexStart = partName[i].Index("WH") + 2;
				sampleProdType = sampleType::WH;
				if(curEvent == 0)
				{
					cout<<"---------------------------"<<endl;
					cout<<"This is a WH sample"<<endl;
					cout<<"---------------------------"<<endl;
				}
			}
			else if(partName[i].Contains("ZH"))
			{
				//cout<<"ZH Sample ";
				indexStart = partName[i].Index("ZH") + 2;
				sampleProdType = sampleType::ZH;
				if(curEvent == 0)
				{
					cout<<"---------------------------"<<endl;
					cout<<"This is a ZH sample"<<endl;
					cout<<"---------------------------"<<endl;
				}
			}
			else if(partName[i].Contains("ttH"))
			{
				//cout<<"ttH Sample ";
				indexStart = partName[i].Index("ttH") + 3;
				sampleProdType = sampleType::ttH;
				if(curEvent == 0)
				{
					cout<<"---------------------------"<<endl;
					cout<<"This is a ttH sample"<<endl;
					cout<<"---------------------------"<<endl;
				}

			}
			else if(partName[i].Contains("qqH"))
			{
				//cout<<"ttH Sample ";
				indexStart = partName[i].Index("qqH") + 3;
				sampleProdType = sampleType::qqF;
				if(curEvent == 0)
				{
					cout<<"---------------------------"<<endl;
					cout<<"This is a qqF sample"<<endl;
					cout<<"---------------------------"<<endl;
				}

			}
			else
			{
				//cout<< "Not Higgs ";
				sampleProdType = sampleType::Background;
				if(curEvent == 0)
				{
					cout<<"---------------------------"<<endl;
					cout<<"This is a not Higgs, so a Background sample"<<endl;
					cout<<"---------------------------"<<endl;
				}
				partName.clear();
				delete parts;
				return -1;
			}
			// Now looping over the name to find where the mass term ends 
			TString nameSample = partName[i];

			Int_t indexEnd = indexStart;
			Int_t len = 0;
			TString charac = nameSample[indexEnd];
			do
			{
				indexEnd++;
				len++;
				charac = nameSample[indexEnd];
			}while(charac.IsDec());
			
			// For the decimal place
			Double_t decimalMass = 0;

			if(nameSample[indexEnd] == 'p')
			{
				//cout<<"There is a point"<<endl;
				Int_t indexEndDec = indexEnd + 1;
				Double_t lenDec = 0;
				TString characDec = nameSample[indexEndDec];
				do
				{
					indexEndDec++;
					lenDec++;
					characDec = nameSample[indexEndDec];
				}while(characDec.IsDec());
				TString higgsMassDecStr = nameSample(indexEnd + 1, lenDec);
				decimalMass = higgsMassDecStr.Atoi();

				// Fix the decimal place
				decimalMass = decimalMass/pow(10.0, lenDec);
				//cout<<"Decimal Mass Higgs: "<<decimalMass<<endl;
			}

			// Getting the mass and returning it
			TString higgsMassStr = nameSample(indexStart, len);
			Int_t higgsMass = higgsMassStr.Atoi();

			Double_t finalHiggsMass = higgsMass+decimalMass;
			if(curEvent == 0)
			{
				cout<<"---------------------------"<<endl;
				cout<<"Sample Mass: "<< finalHiggsMass<<endl;
				cout<<"---------------------------"<<endl;
			}
			partName.clear();
			delete parts;

			return finalHiggsMass;
		}
	}
	// Default return
	// For Debugging purposes
	return -2;
}

void HiggsAnalysis::FillProductionTag()
{
	// Storing the name of the file
	TString name;
	name = currFileName;
	// Overwrite from grid
	if(runningGrid) name = gridFileName;

	// data doesn't have this tag, so just return generic ones
	if(!isMC)
	{
		// if data
		curMCCollection = -1;

		if(dataYear == 2011 && !useNewGeoData) curDataCollection = dataCalibType::y2011c;
		else if(dataYear == 2012 && !useNewGeoData) curDataCollection = dataCalibType::y2012ab;
		else if(dataYear == 2011 && useNewGeoData) curDataCollection = dataCalibType::y2011d;
		else if(dataYear == 2012 && useNewGeoData) curDataCollection = dataCalibType::y2012c;

		return;
	}
	
	// If sample is Z(p)Zp, special case, as it doesn't have s tag
	if(sampleProdType == sampleType::ggF_ZpZp)
	{
		if(dataYear == 2011) curMCCollection = MCCollection::MC11c;
		else if(dataYear == 2012)
		{
			if(event->eventinfo.RunNumber() == 195847) {curMCCollection = MCCollection::MC12a;}
			else if(event->eventinfo.RunNumber() == 195848) {curMCCollection = MCCollection::MC12b;}
			else{cout<<"Error: FillProductionTag: Mc12ab runnumber not recognized"<<endl;}
		}

		return;
	}

	// If not data
	curDataCollection = -1;
	// Vector listing the tags for the production
	Int_t mc11c[] 	= {1272,1273,1274,1299,1300,1309,1310,1349,1350,1351,1352,1353,1370,1372,1378,1571};
	Int_t mc11d[] 	= {1786};
	Int_t mc12ab[] 	= {1468,1469,1470,1472,1479,1482,1484,1485,1486,1499,1504,1581,1586,1589,1599,1609,1610,1611,1716,1773,1773,1776};	
	Int_t mc12c[] 	= {1737,1741,1746,1748,1771,1798,1799,1831,1832};

	// Splitting the file path
	TObjArray *parts = name.Tokenize(".");
	vector<TString> partName;
	if(parts->GetEntriesFast()) {
	   TIter iString(parts);
	   TObjString* os=0;
	   while ((os=(TObjString*)iString())) {
	  	partName.push_back(os->GetString());
	   }
	}

	for(Int_t i = 0; i < (Int_t) partName.size(); i++)
	{
		if(partName[i].Contains("_s"))
		{
			// Splitting part based on "_"
			TObjArray *parts = partName[i].Tokenize("_");
			vector<TString> prodTag;
			if(parts->GetEntriesFast()) {
	   			TIter iString(parts);
	   			TObjString* os=0;
	   			while ((os=(TObjString*)iString())) {
	  			prodTag.push_back(os->GetString());
	   			}
			}
			
			// Production tag
			TString sProdTag = prodTag[1];
			// Just getting the number part
			sProdTag = sProdTag(1,4);
			Int_t prodTagNum = sProdTag.Atoi();

			// Checking m11c
			Int_t nMC11c  =	(sizeof(mc11c)/sizeof(*mc11c));
			Int_t nMC11d  =	(sizeof(mc11d)/sizeof(*mc11d));
			Int_t nMC12ab =	(sizeof(mc12ab)/sizeof(*mc12ab));
			Int_t nMC12c  =	(sizeof(mc12c)/sizeof(*mc12c));

			// Checking what sample it in
			if(find(mc11c,mc11c+nMC11c,prodTagNum) != mc11c+nMC11c)
			{
				curMCCollection = MCCollection::MC11c;
			}
			else if(find(mc11d,mc11d+nMC11d,prodTagNum) != mc11d+nMC11d)
			{
				curMCCollection = MCCollection::MC11d;
			}
			else if(find(mc12ab,mc12ab+nMC12ab,prodTagNum) != mc12ab+nMC12ab)
			{
				if(event->eventinfo.RunNumber() == 195847) {curMCCollection = MCCollection::MC12a;}
				else if(event->eventinfo.RunNumber() == 195848) {curMCCollection = MCCollection::MC12b;}
				else{cout<<"Error: FillProductionTag: Mc12ab runnumber not recognized"<<endl;}
			}
			else if(find(mc12c,mc12c+nMC12c,prodTagNum) != mc12c+nMC12c)
			{
				curMCCollection = MCCollection::MC12c;
			}
			else
			{
				cout<<"Error: FillProductionTag: production tag not recognized"<<endl;
			}

			return;
		}
	}
}

// Gets the trigger name in the form of a strng
void HiggsAnalysis::FillTriggerString(TString singleMu [], TString diMu[], TString singleEl [], TString diEl [], TString eMu[])
{
	// Cleaning the Variables
	for(Int_t i = 0; i < 2; i++)
	{
		singleMu[i] = "";
		diMu[i] = "";
		singleEl[i] = "";
		diEl[i] = "";
		eMu[i] = "";
	}
	// Filling the variables based on Data Period
	if(dataPeriod == period2011_BD || dataPeriod == period2011_EH)
	{
		singleMu[0] = "EF_mu18_MG";
		diMu[0] = "EF_2mu10_loose";

		singleEl[0] = "EF_e20_medium";
		diEl[0] = "EF_2e12_medium";

		eMu[0] = "EF_e10_medium_mu6";
	}
	else if(dataPeriod == period2011_IK)
	{
		// Single Muon Trigger
		if(runNumber_sf < 186516) { singleMu[0] = "EF_mu18_MG"; }
		else { singleMu[0] = "EF_mu18_MG_medium"; }

		// Single and Di Electron Trigger
		if(runNumber_sf < 186873)
		{
			singleEl[0] = "EF_e20_medium";
			diEl[0] = "EF_2e12_medium";
		}
		else
		{
			singleEl[0] = "EF_e22_medium";
			diEl[0] = "EF_2e12T_medium";
		}
		// Other triggers
		diMu[0] = "EF_2mu10_loose";
		eMu[0] = "EF_e10_medium_mu6";
	}
	else if(dataPeriod == period2011_LM)
	{
		singleMu[0] = "EF_mu18_MG_medium";
		diMu[0] = "EF_2mu10_loose";

		singleEl[0] = "EF_e22vh_medium1";
		diEl[0] = "EF_2e12Tvh_medium";

		eMu[0] = "EF_e10_medium_mu6";
	}
	else if(dataPeriod == period2012_All)
	{
		singleMu[0] = "EF_mu24i_tight";
		singleMu[1] = "EF_mu36_tight";
		
		diMu[0] = "EF_2mu13";
		diMu[1] = "EF_mu18_tight_mu8_EFFS";

		singleEl[0] = "EF_e24vhi_medium1";
		singleEl[1] = "EF_e60_medium1";
		
		diEl[0] = "EF_2e12Tvh_loose1";
		if(!isMC) diEl[1] = "EF_2e12Tvh_loose1_L2StarB";

		eMu[0] = "EF_e12Tvh_medium1_mu8";
		eMu[1] = "EF_e24vhi_loose1_mu8";
	}
}
////////////////////////////////////////////////////////////////////////////////////////
//				Printing and Saving...
////////////////////////////////////////////////////////////////////////////////////////
// To print the initial vars
void HiggsAnalysis::PrintInitVar()
{
	cout<<"Intial Variable for the cutFlow Analysis"<<endl;
	cout<<"Data file Name: "<<physicsTree->GetDirectory()->GetName()<<endl;
	cout<<"Run Year: "<< dataYear<<endl;
	cout<<"CM: "<<CM_E<<endl;
	cout<<"is MC: "<<isMC<<endl;
	cout<<"period: "<<dataPeriod<<endl;
 
}

// Debug Functions
// Prints the informations for the given Event Number
void HiggsAnalysis::printDebug(int EntNumber, Bool_t printllee)
{
	if (printllee) isDebugCallllee = true;
	else isDebugCall = true;
	//printWeight = true;
	// For the main Loop
	Int_t nEvent = physicsTree->GetEntries();
	// Main loop
	for(Long64_t iEvent = 0; iEvent < nEvent; iEvent++)
	{
		curEvent = iEvent;
		Long64_t currEvent = iEvent;
		TChain* chain = dynamic_cast<TChain *> (physicsTree);
		if(chain)
		{
			currEvent = chain->LoadTree(currEvent);
			curEvent = currEvent;
		}
		event->GetEntry(currEvent);
		if((int)event->eventinfo.EventNumber() == EntNumber){
			cout<<endl<<endl<<endl;
		cout<<"Print Debug"<<endl;
		cout<<"------------------------"<<endl;
		cout<<"------------------------"<<endl;
		cout<<"------------------------"<<endl;
		
		printDebugInfo();
		return;

		}
	}
}
// Prints the information for this event
void HiggsAnalysis::printDebugInfo()
{
	cout<<"EventNumber: "<<event->eventinfo.EventNumber()<<endl;
	cout<<"RunNumber: "<<event->eventinfo.RunNumber()<<endl;
	getPeriodEvent();
	cout<<"DataPeriod: "<<dataPeriod<<endl;
	cout<<"isMC: "<<isMC<<endl;
	cout<<"Year: "<<dataYear<<endl;
	cout<<"AnalyseEventOutput: (don't really trust this..): "<< AnalyzeTreeEvent(curEvent)<<endl;
	cout<<"DataPreselectionCut: "<<DataPreselectionCut()<<endl;
	cout<<"AllPreselectionCut: "<<AllPreselectionCut()<<endl;
	cout<<"Vertex Cut: "<< VertexCut()<<endl;
	cout<<"Single Electron Trigger: "<<SingleElectronTrigger(runNumber_sf)<<endl;
	cout<<"Di Electron Trigger: "<<DiElectronTrigger(runNumber_sf)<<endl;
	cout<<"Single Muon Trigger: "<<SingleMuonTrigger(runNumber_sf)<<endl;
	cout<<"Di Muon Trigger: "<<DiMuonTrigger()<<endl;
	cout<<"Electron Muon Trigger: "<<ElectronMuonTrigger()<<endl;

	// Fill the counting Hist
	FillCountingHist();
}
// Print MC information
void HiggsAnalysis::printMCInfo(int EntNumber)
{
	// For the main Loop
	Int_t nEvent = physicsTree->GetEntries();
	// Main loop
	for(Long64_t iEvent = 0; iEvent < nEvent; iEvent++)
	{
		curEvent = iEvent;
		Long64_t currEvent = iEvent;
		TChain* chain = dynamic_cast<TChain *> (physicsTree);
		if(chain)
		{
			currEvent = chain->LoadTree(currEvent);
			curEvent = currEvent;
		}
		event->GetEntry(currEvent);
		if((int)event->eventinfo.EventNumber() != EntNumber) continue;
		cout<<"MC information for eventNumber: "<<EntNumber<<endl;
		cout<<"Index\tparticleName\tPDG ID\tMcStatus\tbarcode\tparentSize\tparentIndex\tchildSize\tchildIndex\tPt\teta\tphi\tm"<<endl;
		// Printing the things
		for(Int_t i = 0; i < event->mc.n(); i++)
		{
			cout<<i<<"\t";
			cout<<getParticleName(event->mc[i].pdgId())<<"\t";
			cout<<event->mc[i].pdgId()<<"\t";
			cout<<event->mc[i].status()<<"\t";
			cout<<event->mc[i].barcode()<<"\t";			
			cout<<event->mc[i].parent_index().size()<<"\t";			
			if(event->mc[i].parent_index().size() > 0) 
			{
				for(Int_t j = 0; j < (Int_t) event->mc[i].parent_index().size(); j++) cout<<event->mc[i].parent_index().at(j)<<" ";
				cout<<"\t";
			}
			else cout<<"0\t";

			cout<<event->mc[i].child_index().size()<<"\t";			
			if(event->mc[i].child_index().size() > 0) 
			{
				for(Int_t j = 0; j < (Int_t) event->mc[i].child_index().size(); j++) cout<<event->mc[i].child_index().at(j)<<" ";
				cout<<"\t";
			}
			else cout<<"0\t";
			cout<<event->mc[i].pt()<<"\t";
			cout<<event->mc[i].eta()<<"\t";
			cout<<event->mc[i].phi()<<"\t";
			cout<<event->mc[i].m()<<"\t";			
			cout<<endl;
		}

		// Print the turth muon info
		cout<<endl<<endl<<endl<<"Truth and reco level Info of the muons staco"<<endl;
		cout<<"Index\ttype\ttypeText\tbarCode\tdr\tisMatched\tmothertype\tmothertypeText\tmotherBarcode\tpt\trecoPt\tTrutheta\tTruthphi"<<endl;
		for(Int_t i = 0; i < event->mu_staco.n(); i++)
		{
			cout<<i<<"\t";
			cout<<event->mu_staco[i].truth_type()<<"\t";
			cout<<getParticleName(event->mu_staco[i].truth_type())<<"\t";
			cout<<event->mu_staco[i].truth_barcode()<<"\t";			
			cout<<event->mu_staco[i].truth_dr()<<"\t";
			cout<<event->mu_staco[i].truth_matched()<<"\t";
			cout<<event->mu_staco[i].truth_mothertype()<<"\t";
			cout<<getParticleName(event->mu_staco[i].truth_mothertype())<<"\t";
			cout<<event->mu_staco[i].truth_motherbarcode()<<"\t";			
			cout<<event->mu_staco[i].truth_pt()<<"\t";
			cout<<event->mu_staco[i].pt()<<"\t";
			cout<<event->mu_staco[i].truth_eta()<<"\t";			
			cout<<event->mu_staco[i].truth_phi()<<"\t";			
			cout<<endl;
		}
		// Print the turth muon info
		cout<<endl<<endl<<endl<<"Truth and reco level Info of the muons calo"<<endl;
		cout<<"Index\ttype\ttypeText\tbarCode\tdr\tisMatched\tmothertype\tmothertypeText\tmotherBarcode\tpt\trecoPt\tTrutheta\tTruthphi"<<endl;
		for(Int_t i = 0; i < event->mu_calo.n(); i++)
		{
			cout<<i<<"\t";
			cout<<event->mu_calo[i].truth_type()<<"\t";
			cout<<getParticleName(event->mu_calo[i].truth_type())<<"\t";
			cout<<event->mu_calo[i].truth_barcode()<<"\t";			
			cout<<event->mu_calo[i].truth_dr()<<"\t";
			cout<<event->mu_calo[i].truth_matched()<<"\t";
			cout<<event->mu_calo[i].truth_mothertype()<<"\t";
			cout<<getParticleName(event->mu_calo[i].truth_mothertype())<<"\t";
			cout<<event->mu_calo[i].truth_motherbarcode()<<"\t";			
			cout<<event->mu_calo[i].truth_pt()<<"\t";
			cout<<event->mu_calo[i].pt()<<"\t";
			cout<<event->mu_calo[i].truth_eta()<<"\t";			
			cout<<event->mu_calo[i].truth_phi()<<"\t";			
			cout<<endl;
		}

		// Print the turth electron info
		cout<<endl<<endl<<endl<<"Truth and reco level Info of the electron"<<endl;
		cout<<"Index\ttype\ttypeText\tbarcode\tisMatched\tmothertype\tmothertypeText\tmotherBarcode\tindex\tpt\trecoPt\tTrutheta\tTruthphi"<<endl;
		for(Int_t i = 0; i < el_cur->n(); i++)
		{
			cout<<i<<"\t";
			cout<<(&((*el_cur)[i]))->truth_type()<<"\t";
			cout<<getParticleName((&((*el_cur)[i]))->truth_type())<<"\t";
			cout<<(&((*el_cur)[i]))->truth_barcode()<<"\t";			
			cout<<(&((*el_cur)[i]))->truth_matched()<<"\t";
			cout<<(&((*el_cur)[i]))->truth_mothertype()<<"\t";
			cout<<getParticleName((&((*el_cur)[i]))->truth_mothertype())<<"\t";
			cout<<(&((*el_cur)[i]))->truth_motherbarcode()<<"\t";			
			cout<<(&((*el_cur)[i]))->truth_index()<<"\t";			
			cout<<(&((*el_cur)[i]))->truth_pt()<<"\t";
			cout<<(&((*el_cur)[i]))->pt()<<"\t";
			cout<<(&((*el_cur)[i]))->truth_eta()<<"\t";			
			cout<<(&((*el_cur)[i]))->truth_phi()<<"\t";
			cout<<endl;
		}

		// Print the turth photon info
		cout<<endl<<endl<<endl<<"Truth and reco level Info of the photon"<<endl;
		cout<<"Index\ttype\ttypeText\tbarcode\tisMatched\tmothertype\tmothertypeText\tmotherbarcode\tindex\tpt\trecoPt\tTrutheta\tTruthphi"<<endl;
		for(Int_t i = 0; i < event->ph.n(); i++)
		{
			cout<<i<<"\t";
			cout<<event->ph[i].truth_type()<<"\t";
			cout<<getParticleName(event->ph[i].truth_type())<<"\t";
			cout<<event->ph[i].truth_barcode()<<"\t";			
			cout<<event->ph[i].truth_matched()<<"\t";
			cout<<event->ph[i].truth_mothertype()<<"\t";
			cout<<getParticleName(event->ph[i].truth_mothertype())<<"\t";
			cout<<event->ph[i].truth_motherbarcode()<<"\t";
			cout<<event->ph[i].truth_index()<<"\t";
			cout<<event->ph[i].truth_pt()<<"\t";			
			cout<<event->ph[i].pt()<<"\t";	
			cout<<event->ph[i].truth_eta()<<"\t";			
			cout<<event->ph[i].truth_phi()<<"\t";		
			cout<<endl;
		}
		return;
	}
}
void HiggsAnalysis::printMCEventInfo(QuadLepton * higgs)
{
	//Finding the 2 Z's
	vector<Int_t> zIndex;
	for(Int_t i = 0; i < event->mc.n(); i++)
	{
		if(event->mc[i].pdgId() == 23) zIndex.push_back(i);
	}

	// Sanity check
	if(zIndex.size() != 2) cout<<"ERROR: "<<zIndex.size()<<" ZBoson found"<<endl;

	// Check if the parent of the z is a higgs. Then check is the status 62
	// Store the index of the parent of the z as well
	vector<Int_t> higgsIndex;
	for(Int_t i = 0; i < (Int_t) zIndex.size(); i++)
	{
		Int_t currI = zIndex[i];
		if(event->mc[currI].parent_index().size() <= 0){ cout<<"ERROR: No parent for the Zboson"<<endl; continue;}

		for(Int_t j = 0; j < (Int_t) event->mc[currI].parent_index().size(); j++)
		{
			Int_t parentIndex = event->mc[currI].parent_index().at(j);
			if(event->mc[parentIndex].pdgId() != 25) cout<<"ERROR: Parent of Z is not a higgs"<<endl;
			if(event->mc[parentIndex].status() != 62) cout<<"WARNING: Parent Higgs doesn't have status 62"<<endl;
			//if(event->mc[parentIndex].pdgId() == 25 && event->mc[parentIndex].status() == 62) cout<<"Z parent is a 'Normal' higgs"<<endl;
			higgsIndex.push_back(parentIndex);
		}
	}
	// Getting the unique in vector	
	higgsIndex.erase(std::unique(higgsIndex.begin(), higgsIndex.end()), higgsIndex.end());
	// Sanity checl
	if(higgsIndex.size() != 1) cout<<"ERROR: "<<higgsIndex.size()<<" parent higgs found"<<endl;

	// Getting the child index
	vector<Int_t> leptonIndex;
	for(Int_t i = 0; i < (Int_t) zIndex.size(); i++)
	{
		Int_t currI = zIndex[i];
		if(event->mc[currI].child_index().size() <= 0){ cout<<"ERROR: No child for the Zboson"<<endl; continue;}

		for(Int_t j = 0; j < (Int_t) event->mc[currI].child_index().size(); j++)
		{
			Int_t childIndex = event->mc[currI].child_index().at(j);
			leptonIndex.push_back(childIndex);
		}
	}

	// Getting the leptonChild index
	vector<Int_t> leptonChildIndex;
	for(Int_t i = 0; i < (Int_t) leptonIndex.size(); i++)
	{
		Int_t currI = leptonIndex[i];
		for(Int_t j = 0; j < (Int_t) event->mc[currI].child_index().size(); j++)
		{
			Int_t childIndex = event->mc[currI].child_index().at(j);
			leptonChildIndex.push_back(childIndex);
		}
	}

	
	// Printing the information
	cout<<"EventNumber: "<<event->eventinfo.EventNumber()<<endl;
	cout<<"Index\tparticleName\tPDG ID\tMcStatus\tbarcode\tparentSize\tparentIndex\tchildSize\tchildIndex\tPt\teta\tphi\tm"<<endl;
	for(Int_t i = 0; i < (Int_t) higgsIndex.size(); i++){printMCParticleInfo(higgsIndex[i]);}
	for(Int_t i = 0; i < (Int_t) zIndex.size(); i++){printMCParticleInfo(zIndex[i]);}
	for(Int_t i = 0; i < (Int_t) leptonIndex.size(); i++){printMCParticleInfo(leptonIndex[i]);}
	for(Int_t i = 0; i < (Int_t) leptonChildIndex.size(); i++){printMCParticleInfo(leptonChildIndex[i]);}

	vector<ChargedLepton *> lep = higgs->getLepton();
	for(vector<ChargedLepton *>::iterator itr_lep = lep.begin();
		itr_lep != lep.end(); ++itr_lep)
	{
		ChargedLepton* curr_i = *itr_lep;
		if(curr_i->getFlavor() == flavor::Electron)
		{
			D3PDReader::ElectronD3PDObjectElement* el =  curr_i->GetElectron();
			cout<<"Electron type\ttypeText\tbarcode\tmothertype\tmothertypeText\tmotherBarcode\tindex\tpt\tTrutheta\tTruthphi"<<endl;
		
			cout<<el->truth_type()<<"\t";
			cout<<getParticleName(el->truth_type())<<"\t";
			cout<<el->truth_barcode()<<"\t";			
			cout<<el->truth_mothertype()<<"\t";
			cout<<getParticleName(el->truth_mothertype())<<"\t";
			cout<<el->truth_motherbarcode()<<"\t";			
			cout<<el->truth_index()<<"\t";			
			cout<<el->truth_pt()<<"\t";
			cout<<el->truth_eta()<<"\t";			
			cout<<el->truth_phi()<<"\t";
			cout<<endl;

		}
		else if(curr_i->getFlavor() == flavor::Muon)
		{
			D3PDReader::MuonD3PDObjectElement* mu =  curr_i->GetMuon();

			cout<<"Muon type\ttypeText\tbarCode\tmothertype\tmothertypeText\tmotherBarcode\tpt\tTrutheta\tTruthphi"<<endl;
			cout<<mu->truth_type()<<"\t";
			cout<<getParticleName(mu->truth_type())<<"\t";
			cout<<mu->truth_barcode()<<"\t";			
			cout<<mu->truth_mothertype()<<"\t";
			cout<<getParticleName(mu->truth_mothertype())<<"\t";
			cout<<mu->truth_motherbarcode()<<"\t";			
			cout<<mu->truth_pt()<<"\t";
			cout<<mu->truth_eta()<<"\t";			
			cout<<mu->truth_phi()<<"\t";			
			cout<<endl;
		}
	}

	cout<<endl;
	
	

}
void HiggsAnalysis::printMCParticleInfo(Int_t index)
{
		// Printing the things
	Int_t i = index;
	cout<<i<<"\t";
	cout<<getParticleName(event->mc[i].pdgId())<<"\t";
	cout<<event->mc[i].pdgId()<<"\t";
	cout<<event->mc[i].status()<<"\t";
	cout<<event->mc[i].barcode()<<"\t";			
	cout<<event->mc[i].parent_index().size()<<"\t";			
	if(event->mc[i].parent_index().size() > 0) 
	{
		for(Int_t j = 0; j < (Int_t) event->mc[i].parent_index().size(); j++) cout<<event->mc[i].parent_index().at(j)<<" ";
		cout<<"\t";
	}
	else cout<<"0\t";

	cout<<event->mc[i].child_index().size()<<"\t";			
	if(event->mc[i].child_index().size() > 0) 
	{
		for(Int_t j = 0; j < (Int_t) event->mc[i].child_index().size(); j++) cout<<event->mc[i].child_index().at(j)<<" ";
		cout<<"\t";
	}
	else cout<<"0\t";
	cout<<event->mc[i].pt()<<"\t";
	cout<<event->mc[i].eta()<<"\t";
	cout<<event->mc[i].phi()<<"\t";
	cout<<event->mc[i].m()<<"\t";			
	cout<<endl;
}

TString HiggsAnalysis::getParticleName(int pdgID)
{
	pdgID = TMath::Abs(pdgID);
	switch (pdgID)
	{
		case 1: return "d";
		case 2: return "u";
		case 3: return "s";
		case 4: return "c";
		case 5: return "b";
		case 6: return "t";
		
		case 11: return "El";
		case 12: return "ElNu";
		case 13: return "Mu";
		case 14: return "MuNu";
		case 15: return "Tau";
		case 16: return "TauNu";

		case 21: return "gluon";
		case 22: return "photon";
		case 23: return "Z_boson";
		case 24: return "w_boson";
		case 25: return "Higgs";

		case 211: return "pion";
		case 111: return "pion";
		case 130: return "kaon";
		
		default: return "Other";
	}
}
////////////////////////////////////////////////////////////////////////////////////////
//				Initialize the vars...
////////////////////////////////////////////////////////////////////////////////////////

// Initialize all the variables necessary for the cutflow
int HiggsAnalysis::InitializeVar()
{
		// Cut Flow the Overall Analysis	
	// For Counting the Event
	nCut = 6;
	cutName = new TString [nCut];
	cutPass = new Int_t [nCut];
	cutPassW = new Double_t [nCut];
	
	cutName[cutFlow::Total] = "Total";
	cutName[cutFlow::DataPreselection] = "DataPreSelection";
	cutName[cutFlow::Preselection] = "PreSelection";
	cutName[cutFlow::Trigger] = "Trigger Main";
	cutName[cutFlow::Trigger4Mu] = "Trigger 4Mu";
	cutName[cutFlow::Trigger4e] = "Trigger 4e";
	for(Int_t i = 0; i < nCut; i++) {cutPass[i] = 0; cutPassW[i] = 0;}
	
	// Cut Flow for the Muon
	nMuCut = 14;
	cutMuName = new TString [nMuCut];
	cutMuPass = new Int_t [nMuCut];

	cutMuName[cutMuFlow::Total] = "Total";
	cutMuName[cutMuFlow::DataPreselection] = "DataPreSelection";
	cutMuName[cutMuFlow::Preselection] = "PreSelection";
	cutMuName[cutMuFlow::Trigger] = "Trigger";
	cutMuName[cutMuFlow::Author] = "Author";
	cutMuName[cutMuFlow::Pt] = "Pt";
	cutMuName[cutMuFlow::Eta] = "Eta";
	cutMuName[cutMuFlow::BLayer] = "BLayer";
	cutMuName[cutMuFlow::Pix] = "Pix";
	cutMuName[cutMuFlow::SCT] = "SCT";
	cutMuName[cutMuFlow::Holes] = "Holes";
	cutMuName[cutMuFlow::TRT] = "TRT";
	cutMuName[cutMuFlow::D0] = "D0/Z0";
	cutMuName[cutMuFlow::OverLap] = "OverLap";
	for(Int_t i = 0; i < nMuCut; i++) {cutMuPass[i] = 0;}

	// Cut Flow for the Electron 
	nElCut = 13;
	cutElName = new TString [nElCut];
	cutElPass = new Int_t [nElCut];
	cutElLoosePass = new Int_t [nElCut];

	cutElName[cutElFlow::Total] = "Total";
	cutElName[cutElFlow::DataPreselection] = "DataPreSelection";
	cutElName[cutElFlow::Preselection] = "PreSelection";
	cutElName[cutElFlow::Trigger] = "Trigger";
	cutElName[cutElFlow::Author] = "Author";
	cutElName[cutElFlow::Loose] = "Loose";
	cutElName[cutElFlow::Eta] = "Eta";
	cutElName[cutElFlow::Et] = "Et";
	cutElName[cutElFlow::ObjectQuality] = "Object Quality";
	cutElName[cutElFlow::Z0] = "Z0";
	cutElName[cutElFlow::OverLapElEl] = "OverLapElEl";
	cutElName[cutElFlow::OverLapClElEl] = "OverLapClElEl";
	cutElName[cutElFlow::OverLap] = "OverLap";
	for(Int_t i = 0; i < nElCut; i++) {cutElPass[i] = 0;}
	for(Int_t i = 0; i < nElCut; i++) {cutElLoosePass[i] = 0;}
	
	// Cut Flow for jets
	nJetsCut = 9;
	cutJetsName = new TString [nJetsCut];
	cutJetsPass = new Int_t [nJetsCut];

	cutJetsName[cutJetsFlow::Total] = "Total";
	cutJetsName[cutJetsFlow::DataPreselection] = "DataPreSelection";
	cutJetsName[cutJetsFlow::Preselection] = "PreSelection";
	cutJetsName[cutJetsFlow::Trigger] = "Trigger";
	cutJetsName[cutJetsFlow::Pt] = "Pt";
	cutJetsName[cutJetsFlow::Eta] = "Eta";
	cutJetsName[cutJetsFlow::Pileup] = "Pileup";
	cutJetsName[cutJetsFlow::Clean] = "Clean";
	cutJetsName[cutJetsFlow::OverLap] = "OverLap";
	
	for(Int_t i = 0; i < nJetsCut; i++) {cutJetsPass[i] = 0;}

	// CutFlow for each channel
	nCH = 13;
	cutCHName = new TString [nCH];
	cut4MuPass = new Int_t [nCH];
	cut4ElPass = new Int_t [nCH];
	cut2L2LPass = new Int_t [nCH];
	cutlleePass = new Int_t [nCH];
	cut4MuPassW = new Double_t [nCH];
	cut4ElPassW = new Double_t [nCH];
	cut2L2LPassW = new Double_t [nCH];
	
	cutCHName[cutFlowCH::Total] = "Total";
	cutCHName[cutFlowCH::Trigger] = "Trigger";
	cutCHName[cutFlowCH::Lepton] = "Lepton";	
	cutCHName[cutFlowCH::SFOS] = "SFOS";
	cutCHName[cutFlowCH::Kinematics] = "Kinematics";
	cutCHName[cutFlowCH::TriggerMatch] = "TriggerMatch";
	cutCHName[cutFlowCH::Z1Mass] = "Z1Mass";
	cutCHName[cutFlowCH::Z2Mass] = "Z2Mass";
	cutCHName[cutFlowCH::DeltaR] = "DeltaR";
	cutCHName[cutFlowCH::TrackIso] = "TrackIso";
	cutCHName[cutFlowCH::CaloIso] = "CaloIso";
	cutCHName[cutFlowCH::D0Sig] = "D0Sig";
	cutCHName[cutFlowCH::Final] = "Final";

	for(Int_t i = 0; i < nCH; i++)
	{
		cut4MuPass[i] = 0; cut4MuPassW[i] = 0;
		cut4ElPass[i] = 0; cut4ElPassW[i] = 0;
		cut2L2LPass[i] = 0; cut2L2LPassW[i] = 0;
		cutlleePass[i] = 0;
	}

	nProdCH = 5;
	prodCHName = new TString [nProdCH];
	prodCH4Mu = new Int_t [nProdCH];
	prodCH4El = new Int_t [nProdCH];
	prodCH2L2L = new Int_t [nProdCH];
	prodCHllee = new Int_t [nProdCH];

	prodCHName [productionChannel::VBF] = "VBF";
	prodCHName [productionChannel::VHLep] = "VHLep";
	prodCHName [productionChannel::VHHad] = "VHHad";
	prodCHName [productionChannel::ggF] = "ggF";
	prodCHName [productionChannel::VH] = "VH";
	for(Int_t i = 0; i < nProdCH; i++)
	{
		prodCH4Mu[i] = 0;
		prodCH4El[i] = 0; 
		prodCH2L2L[i] = 0;
		prodCHllee[i] = 0;
	}
	
	nListTruthQuadType = 7;
	truthQuadType = new TString [nListTruthQuadType];
	nTruthQuadType = new Int_t [nListTruthQuadType];
	for(Int_t i = 0; i < nListTruthQuadType; i++) {nTruthQuadType[i] = 0;}
	truthQuadType[truthTypeQuad::_4mu] 		= "4mu";
	truthQuadType[truthTypeQuad::_4e] 		= "4e";
	truthQuadType[truthTypeQuad::_4tau] 	= "4tau";
	truthQuadType[truthTypeQuad::_2mu2e]	= "2mu2e";
	truthQuadType[truthTypeQuad::_2tau2mu] 	= "2tau2mu";
	truthQuadType[truthTypeQuad::_2tau2e]	= "2tau2e";
	truthQuadType[truthTypeQuad::noStatus]	= "noStatus";

		
	return 0;
}

// Initializing the tools for trigger matching
// Essentially copied from Fabien
void HiggsAnalysis::InitTriggerMatchingTool()
{
    triggerNavigationVariables->set_trig_DB_SMK(event->trig_DB.SMK());
    triggerNavigationVariables->set_trig_Nav_n(event->trig_Nav.n());
    triggerNavigationVariables->set_trig_Nav_chain_ChainId(event->trig_Nav.chain_ChainId());
    triggerNavigationVariables->set_trig_Nav_chain_RoIType(event->trig_Nav.chain_RoIType());
    triggerNavigationVariables->set_trig_Nav_chain_RoIIndex(event->trig_Nav.chain_RoIIndex());

    // electron 
    triggerNavigationVariables->set_trig_RoI_EF_e_egammaContainer_egamma_Electrons(event->trig_RoI_EF_e.egammaContainer_egamma_Electrons());
    triggerNavigationVariables->set_trig_RoI_EF_e_egammaContainer_egamma_ElectronsStatus(event->trig_RoI_EF_e.egammaContainer_egamma_ElectronsStatus());
    triggerNavigationVariables->set_trig_EF_el_n(event->trig_EF_el.n());
    triggerNavigationVariables->set_trig_EF_el_eta(event->trig_EF_el.eta());
    triggerNavigationVariables->set_trig_EF_el_phi(event->trig_EF_el.phi());

    // muon 
    triggerNavigationVariables->set_trig_RoI_EF_mu_Muon_ROI(event->trig_RoI_EF_mu.Muon_ROI());
    if(!D3DPTriggerDev){
      triggerNavigationVariables->set_trig_RoI_EF_mu_TrigMuonEFInfoContainer(event->trig_RoI_EF_mu.TrigMuonEFInfoContainer());
      triggerNavigationVariables->set_trig_RoI_EF_mu_TrigMuonEFInfoContainerStatus(event->trig_RoI_EF_mu.TrigMuonEFInfoContainerStatus());
    }else{
      triggerNavigationVariables->set_trig_RoI_EF_mu_TrigMuonEFInfoContainerStatus(event->trig_RoI_EF_mu.TrigMuonEFInfoContainer_eMuonEFInfoStatus());
      triggerNavigationVariables->set_trig_RoI_EF_mu_TrigMuonEFInfoContainer(event->trig_RoI_EF_mu.TrigMuonEFInfoContainer_eMuonEFInfo());
    }
    if(dataYear==2012){
      // for 2012 isolated trigger
      triggerNavigationVariables->set_trig_RoI_EF_mu_TrigMuonEFIsolationContainer(event->trig_RoI_EF_mu.TrigMuonEFIsolationContainer());
      triggerNavigationVariables->set_trig_RoI_EF_mu_TrigMuonEFIsolationContainerStatus(event->trig_RoI_EF_mu.TrigMuonEFIsolationContainerStatus());
	}
    triggerNavigationVariables->set_trig_RoI_L2_mu_CombinedMuonFeature(event->trig_RoI_L2_mu.CombinedMuonFeature());
    triggerNavigationVariables->set_trig_RoI_L2_mu_CombinedMuonFeatureStatus(event->trig_RoI_L2_mu.CombinedMuonFeatureStatus());
    triggerNavigationVariables->set_trig_RoI_L2_mu_MuonFeature(event->trig_RoI_L2_mu.MuonFeature());
    triggerNavigationVariables->set_trig_RoI_L2_mu_Muon_ROI(event->trig_RoI_L2_mu.Muon_ROI());
    triggerNavigationVariables->set_trig_EF_trigmuonef_track_MuonType(event->trig_EF_trigmuonef.track_MuonType());
    triggerNavigationVariables->set_trig_EF_trigmuonef_track_CB_pt(event->trig_EF_trigmuonef.track_CB_pt());
    triggerNavigationVariables->set_trig_EF_trigmuonef_track_CB_eta(event->trig_EF_trigmuonef.track_CB_eta());  
    triggerNavigationVariables->set_trig_EF_trigmuonef_track_CB_phi(event->trig_EF_trigmuonef.track_CB_phi());
    triggerNavigationVariables->set_trig_EF_trigmuonef_track_SA_pt(event->trig_EF_trigmuonef.track_SA_pt());
    triggerNavigationVariables->set_trig_EF_trigmuonef_track_SA_eta(event->trig_EF_trigmuonef.track_SA_eta());
    triggerNavigationVariables->set_trig_EF_trigmuonef_track_SA_phi(event->trig_EF_trigmuonef.track_SA_phi());
    triggerNavigationVariables->set_trig_EF_trigmugirl_track_CB_pt(event->trig_EF_trigmugirl.track_CB_pt());
    triggerNavigationVariables->set_trig_EF_trigmugirl_track_CB_eta(event->trig_EF_trigmugirl.track_CB_eta());
    triggerNavigationVariables->set_trig_EF_trigmugirl_track_CB_phi(event->trig_EF_trigmugirl.track_CB_phi());
    triggerNavigationVariables->set_trig_L2_combmuonfeature_eta(event->trig_L2_combmuonfeature.eta());
    triggerNavigationVariables->set_trig_L2_combmuonfeature_phi(event->trig_L2_combmuonfeature.phi());
    triggerNavigationVariables->set_trig_L2_muonfeature_eta(event->trig_L2_muonfeature.eta());
    triggerNavigationVariables->set_trig_L2_muonfeature_phi(event->trig_L2_muonfeature.phi());
    triggerNavigationVariables->set_trig_L1_mu_eta(event->trig_L1_mu.eta());
    triggerNavigationVariables->set_trig_L1_mu_phi(event->trig_L1_mu.phi());
    triggerNavigationVariables->set_trig_L1_mu_thrName(event->trig_L1_mu.thrName());
    
	triggerNavigationVariables->set_trig_EF_el_Et(event->trig_EF_el.Et());

    
    if (!D3DPTriggerDev  && !triggerNavigationVariables->isValid()) {
      std::cout <<"Trigger Matching Tool: variables not correctly set !"<< std::endl;
    }
	 if (!triggerNavigationVariables->isValid()) {
      std::cout <<"Trigger Matching Tool: variables not correctly set !"<< std::endl;
    }
        
    muonTriggerMatchTool->setDeltaR(0.15);
    electronTriggerMatchTool->setDeltaR(0.15);

}
void HiggsAnalysis::InitTriggerMatchingToolMain()
{
	triggerNavigationVariables= new TriggerNavigationVariables();
    muonTriggerMatchTool = new MuonTriggerMatching(triggerNavigationVariables);
    electronTriggerMatchTool = new ElectronTriggerMatching(triggerNavigationVariables); 
	if(dataYear == 2011)
	{
		leptonSF = new LeptonTriggerSF(2011, 
                                       "../../TrigMuonEfficiency/share", 
                                       "muon_trigger_sf_mc11c.root",
                                       "../../ElectronEfficiencyCorrection/data/",
                                       "rel17p0.v02");
	}
	else if(dataYear == 2012)
	{
		leptonSF = new LeptonTriggerSF(2012, 
                                       "../../TrigMuonEfficiency/share", 
                                       "muon_trigger_sf_2012_AtoL.p1328.root",
                                       "../../ElectronEfficiencyCorrection/data/",
                                       "rel17p2.v07");
	}
}

// To setup the print list. Called from the steering script
void HiggsAnalysis::SetupPrintEventList(Bool_t overWrite, TString fileName)
{
	if(!printEventList) return;

	if(fileName.Contains(".txt"))
	{
		Int_t index = fileName.Index(".txt", 0);
		fileName = fileName(0, index);
	}
	TString mu4File; 
	TString el4File;
	TString el2mu2File;
	TString mu2el2File;
	TString l2l2File;
	TString sampleNameFile;
	// File name based on if running on grid or not
	if(runningGrid)
	{
		mu4File = fileName+"_4Mu.txt";
		el4File = fileName+"_4El.txt";
		el2mu2File = fileName+"_2El2Mu.txt";
		mu2el2File = fileName+"_2Mu2El.txt";
		l2l2File = fileName+"_2El2Mu_2Mu2El.txt";
		sampleNameFile = fileName+"_sampleName.txt";
	}
	else
	{
		mu4File = fileName+TString::Itoa(dataYear,10)+"_4Mu.txt";
		el4File = fileName+TString::Itoa(dataYear,10)+"_4El.txt";
		el2mu2File = fileName+TString::Itoa(dataYear,10)+"_2El2Mu.txt";
		mu2el2File = fileName+TString::Itoa(dataYear,10)+"_2Mu2El.txt";
		l2l2File = fileName+TString::Itoa(dataYear,10)+"_2El2Mu_2Mu2El.txt";
		sampleNameFile = fileName+TString::Itoa(dataYear,10)+"_sampleName.txt";
	}

	// Creating the files
	if(overWrite)
	{
		file4Mu.open("EventList/"+mu4File, ios::out);
		file4El.open("EventList/"+el4File, ios::out);
		file2El2Mu.open("EventList/"+el2mu2File, ios::out);
		file2Mu2El.open("EventList/"+mu2el2File, ios::out);
		file2L2L.open("EventList/"+l2l2File, ios::out);
		fileSampleName.open("EventList/"+sampleNameFile, ios::out);

		file4Mu<<"H->ZZ->4mu Event list\n";
		file4El<<"H->ZZ->4e Event list\n";
		file2El2Mu<<"H->ZZ->2e2mu Event list\n";
		file2Mu2El<<"H->ZZ->2mu2e Event list\n";
		file2L2L<<"H->ZZ->2e2mu+2mu2e Event list\n";
		fileSampleName<<"Event list\n";

		if(printMass)
		{
			file4Mu<<"EventNumber\tMass4l[GeV]\tMassErr4l[GeV]\tMassZ1[GeV]\tMassZ2[GeV]\tMass4lFSR[GeV]\tMassErr4lFSR[GeV]\t";
			file4Mu<<"MassZ1FSR[GeV]\tMassZ2FSR[GeV]\tMass4lZCont\tMassErr4lZCont\tMassZ1ZCont\tMassZ2ZCont\n";
			file4El<<"EventNumber\tMass4l[GeV]\tMassErr4l[GeV]\tMassZ1[GeV]\tMassZ2[GeV]\tMass4lFSR[GeV]\tMassErr4lFSR[GeV]\t";
			file4El<<"MassZ1FSR[GeV]\tMassZ2FSR[GeV]\tMass4lZCont\tMassErr4lZCont\tMassZ1ZCont\tMassZ2ZCont\n";
			file2El2Mu<<"EventNumber\tMass4l[GeV]\tMassErr4l[GeV]\tMassZ1[GeV]\tMassZ2[GeV]\tMass4lFSR[GeV]\tMassErr4lFSR[GeV]\t";
			file2El2Mu<<"MassZ1FSR[GeV]\tMassZ2FSR[GeV]\tMass4lZCont\tMassErr4lZCont\tMassZ1ZCont\tMassZ2ZCont\n";
			file2Mu2El<<"EventNumber\tMass4l[GeV]\tMassErr4l[GeV]\tMassZ1[GeV]\tMassZ2[GeV]\tMass4lFSR[GeV]\tMassErr4lFSR[GeV]\t";
			file2Mu2El<<"MassZ1FSR[GeV]\tMassZ2FSR[GeV]\tMass4lZCont\tMassErr4lZCont\tMassZ1ZCont\tMassZ2ZCont\n";
			file2L2L<<"EventNumber\tMass4l[GeV]\tMassErr4l[GeV]\tMassZ1[GeV]\tMassZ2[GeV]\tMass4lFSR[GeV]\tMassErr4lFSR[GeV]\t";
			file2L2L<<"MassZ1FSR[GeV]\tMassZ2FSR[GeV]\tMass4lZCont\tMassErr4lZCont\tMassZ1ZCont\tMassZ2ZCont\n";	
			fileSampleName<<"EventNumber\tAnalysisType\tSampleName"<<endl;
		}
	}
	else
	{
		file4Mu.open("EventList/"+mu4File,ios::out | ios::app);
		file4El.open("EventList/"+el4File,ios::out | ios::app);
		file2El2Mu.open("EventList/"+el2mu2File,ios::out | ios::app);	
		file2Mu2El.open("EventList/"+mu2el2File,ios::out | ios::app);
		file2L2L.open("EventList/"+l2l2File,ios::out | ios::app);
	}

}
// Actually outputs the print list
void HiggsAnalysis::PrintEventList(Bool_t passCut4Mu, Bool_t passCut4El, Bool_t passCut2L2L)
{
	// Choosing the Bool stream for 2e2mu and 2mu2e channel
	Bool_t is2El2Mu = true;
	if(passCut2L2L)
	{
		if(higgsCandidate2L2L->getQuadType() == quadType::Mu2El2)
			is2El2Mu = false;
		else if(higgsCandidate2L2L->getQuadType() == quadType::El2Mu2)
			is2El2Mu = true;
		else 
		{
			cout<<"PrintEventList(): 2L2L candidate Not recognized"<<endl;
			cout<<"Default to the 2El2Mu outputstream"<<endl;
			is2El2Mu = true;
		}
	}
	// Printing just the events
	if(passCut4Mu) file4Mu<<event->eventinfo.EventNumber();
	if(passCut4El) file4El<<event->eventinfo.EventNumber();
	if(passCut2L2L && is2El2Mu) file2El2Mu<<event->eventinfo.EventNumber();
	if(passCut2L2L && !is2El2Mu) file2Mu2El<<event->eventinfo.EventNumber();
	if(passCut2L2L) file2L2L<<event->eventinfo.EventNumber();
	
	// Printting the mass if the flag is set
	if(passCut4Mu && printMass) 
	{
		file4Mu<<"\t"<<higgsCandidate4Mu->get4Momentum()->M()/1000<<"\t";
		file4Mu<<higgsCandidate4Mu->getMassErr()/1000<<"\t";
		file4Mu<<higgsCandidate4Mu->getZ1()->get4Momentum()->M()/1000<<"\t";
		file4Mu<<higgsCandidate4Mu->getZ2()->get4Momentum()->M()/1000<<"\t";
		file4Mu<<higgsCandidate4Mu->getMassFSR()/1000<<"\t";
		file4Mu<<higgsCandidate4Mu->getMassErrFSR()/1000<<"\t";
		file4Mu<<higgsCandidate4Mu->getZ1MassFSR()/1000<<"\t";
		file4Mu<<higgsCandidate4Mu->getZ2()->get4Momentum()->M()/1000<<"\t";
		file4Mu<<higgsCandidate4Mu->getMassZMassCons()/1000<<"\t";
		file4Mu<<higgsCandidate4Mu->getMassErrZMassCons()/1000<<"\t";
		file4Mu<<higgsCandidate4Mu->getZ1MassZMassCons()/1000<<"\t";
		file4Mu<<higgsCandidate4Mu->getZ2MassZMassCons()/1000;
	}
	if(passCut4El && printMass) 
	{
		file4El<<"\t"<<higgsCandidate4El->get4Momentum()->M()/1000<<"\t";
		file4El<<higgsCandidate4El->getMassErr()/1000<<"\t";		
		file4El<<higgsCandidate4El->getZ1()->get4Momentum()->M()/1000<<"\t";
		file4El<<higgsCandidate4El->getZ2()->get4Momentum()->M()/1000<<"\t";
		file4El<<higgsCandidate4El->getMassFSR()/1000<<"\t";
		file4El<<higgsCandidate4El->getMassErrFSR()/1000<<"\t";
		file4El<<higgsCandidate4El->getZ1MassFSR()/1000<<"\t";
		file4El<<higgsCandidate4El->getZ2()->get4Momentum()->M()/1000<<"\t";
		file4El<<higgsCandidate4El->getMassZMassCons()/1000<<"\t";
		file4El<<higgsCandidate4El->getMassErrZMassCons()/1000<<"\t";
		file4El<<higgsCandidate4El->getZ1MassZMassCons()/1000<<"\t";
		file4El<<higgsCandidate4El->getZ2MassZMassCons()/1000;
	}
	if(passCut2L2L && is2El2Mu && printMass) 
	{
		file2El2Mu<<"\t"<<higgsCandidate2L2L->get4Momentum()->M()/1000<<"\t";
		file2El2Mu<<higgsCandidate2L2L->getMassErr()/1000<<"\t";		
		file2El2Mu<<higgsCandidate2L2L->getZ1()->get4Momentum()->M()/1000<<"\t";
		file2El2Mu<<higgsCandidate2L2L->getZ2()->get4Momentum()->M()/1000<<"\t";
		file2El2Mu<<higgsCandidate2L2L->getMassFSR()/1000<<"\t";
		file2El2Mu<<higgsCandidate2L2L->getMassErrFSR()/1000<<"\t";
		file2El2Mu<<higgsCandidate2L2L->getZ1MassFSR()/1000<<"\t";
		file2El2Mu<<higgsCandidate2L2L->getZ2()->get4Momentum()->M()/1000<<"\t";
		file2El2Mu<<higgsCandidate2L2L->getMassZMassCons()/1000<<"\t";
		file2El2Mu<<higgsCandidate2L2L->getMassErrZMassCons()/1000<<"\t";
		file2El2Mu<<higgsCandidate2L2L->getZ1MassZMassCons()/1000<<"\t";
		file2El2Mu<<higgsCandidate2L2L->getZ2MassZMassCons()/1000;
	}
	if(passCut2L2L && !is2El2Mu && printMass) 
	{
		file2Mu2El<<"\t"<<higgsCandidate2L2L->get4Momentum()->M()/1000<<"\t";
		file2Mu2El<<higgsCandidate2L2L->getMassErr()/1000<<"\t";		
		file2Mu2El<<higgsCandidate2L2L->getZ1()->get4Momentum()->M()/1000<<"\t";
		file2Mu2El<<higgsCandidate2L2L->getZ2()->get4Momentum()->M()/1000<<"\t";
		file2Mu2El<<higgsCandidate2L2L->getMassFSR()/1000<<"\t";
		file2Mu2El<<higgsCandidate2L2L->getMassErrFSR()/1000<<"\t";
		file2Mu2El<<higgsCandidate2L2L->getZ1MassFSR()/1000<<"\t";
		file2Mu2El<<higgsCandidate2L2L->getZ2()->get4Momentum()->M()/1000<<"\t";
		file2Mu2El<<higgsCandidate2L2L->getMassZMassCons()/1000<<"\t";
		file2Mu2El<<higgsCandidate2L2L->getMassErrZMassCons()/1000<<"\t";
		file2Mu2El<<higgsCandidate2L2L->getZ1MassZMassCons()/1000<<"\t";
		file2Mu2El<<higgsCandidate2L2L->getZ2MassZMassCons()/1000;
	}
	if(passCut2L2L  && printMass) 
	{
		file2L2L<<"\t"<<higgsCandidate2L2L->get4Momentum()->M()/1000<<"\t";
		file2L2L<<higgsCandidate2L2L->getMassErr()/1000<<"\t";		
		file2L2L<<higgsCandidate2L2L->getZ1()->get4Momentum()->M()/1000<<"\t";
		file2L2L<<higgsCandidate2L2L->getZ2()->get4Momentum()->M()/1000<<"\t";
		file2L2L<<higgsCandidate2L2L->getMassFSR()/1000<<"\t";
		file2L2L<<higgsCandidate2L2L->getMassErrFSR()/1000<<"\t";
		file2L2L<<higgsCandidate2L2L->getZ1MassFSR()/1000<<"\t";
		file2L2L<<higgsCandidate2L2L->getZ2()->get4Momentum()->M()/1000<<"\t";
		file2L2L<<higgsCandidate2L2L->getMassZMassCons()/1000<<"\t";
		file2L2L<<higgsCandidate2L2L->getMassErrZMassCons()/1000<<"\t";
		file2L2L<<higgsCandidate2L2L->getZ1MassZMassCons()/1000<<"\t";
		file2L2L<<higgsCandidate2L2L->getZ2MassZMassCons()/1000;
	}	


	// Printing the new line charcter
	if(passCut4Mu) file4Mu<<"\n";
	if(passCut4El) file4El<<"\n";
	if(passCut2L2L && is2El2Mu) file2El2Mu<<"\n";
	if(passCut2L2L && !is2El2Mu) file2Mu2El<<"\n";
	if(passCut2L2L ) file2L2L<<"\n";

	// For SampleName File
	if(passCut4Mu) fileSampleName<<event->eventinfo.EventNumber()<<"\t4Mu\t"<<currFileName;
	if(passCut4El) fileSampleName<<event->eventinfo.EventNumber()<<"\t4El\t"<<currFileName;
	if(passCut2L2L && is2El2Mu) fileSampleName<<event->eventinfo.EventNumber()<<"\t2El2Mu\t"<<currFileName;
	if(passCut2L2L && !is2El2Mu) fileSampleName<<event->eventinfo.EventNumber()<<"\t2Mu2El\t"<<currFileName;

	if(passCut4Mu && runningGrid) fileSampleName<<"\t"<<gridFileName;
	if(passCut4El && runningGrid) fileSampleName<<"\t"<<gridFileName;
	if(passCut2L2L && is2El2Mu && runningGrid) fileSampleName<<"\t"<<gridFileName;
	if(passCut2L2L && !is2El2Mu && runningGrid) fileSampleName<<"\t"<<gridFileName;

	if(passCut4Mu) fileSampleName<<endl;
	if(passCut4El) fileSampleName<<endl;
	if(passCut2L2L && is2El2Mu) fileSampleName<<endl;
	if(passCut2L2L && !is2El2Mu) fileSampleName<<endl;

}

// To save the histograms. Called from outside as well
void HiggsAnalysis::SaveHist(Bool_t overWrite)
{

	// Opening the file
	TFile *output;
	if(runningGrid)
	{
		output = new TFile ("Output/EventSummaryPlots.root", "RECREATE");
	}
	else
	{
		if(overWrite) {output = new TFile ("Output/Events"+TString::Itoa(dataYear,10)+".root", "RECREATE");}
		else {output = new TFile ("Output/Events"+TString::Itoa(dataYear,10)+".root", "UPDATE");}
	}

	if (output->IsZombie()) {cout << "Error opening file" << endl; }

	// Saving the Histrograms
	Hist->SaveHist(output, getSampleName(), countingHist);	
}

// To get the sample name of the data set
TString HiggsAnalysis::getSampleName()
{
	// Getting the path of the current TTree.. so that I can save it in the right folder
	TString fileName;
	fileName = currFileName;
	// Overwrite from grid
	if(runningGrid) {return gridFileName;}

	// Spliting the path based on "/" delimeter
	// Code copied from Valerios
	std::vector<TString> fileNamePart;
	TObjArray *chains = fileName.Tokenize("/");
   	if (chains->GetEntriesFast()) {
      TIter iChain(chains);
      TObjString *os = 0;
	      while ((os = (TObjString*)iChain())) {
         	fileNamePart.push_back(os->GetString().Data());
	      } // loop over chains
	} 
	if(printWeight) return "debug call printWeight";
	return fileNamePart[fileNamePart.size() - 2];
}
// To populate the necessary vars
void HiggsAnalysis::FillCountingHist()
{
	// Filing cutPass
	for(Int_t i = 0; i < nCut; i++)
	 {Hist->cutPassHist->SetBinContent(i+1, cutPass[i]); 
	 cutPassW[i] = Hist->cutPassHistW->GetBinContent(i+1);}

	// Filing cutMuPass
	for(Int_t i = 0; i < nMuCut; i++)
	 {Hist->cutMuPassHist->SetBinContent(i+1, cutMuPass[i]); }
	// Filing cutElPass
	for(Int_t i = 0; i < nElCut; i++)
	 {Hist->cutElPassHist->SetBinContent(i+1, cutElPass[i]); }
	// Filing cutJetsPass
	for(Int_t i = 0; i < nJetsCut; i++)
	 {Hist->cutJetsPassHist->SetBinContent(i+1, cutJetsPass[i]);}
	// Filing QuadLepton
	for(Int_t i = 0; i < nCH; i++)
	 {
		Hist->cut4MuPassHist->SetBinContent(i+1, cut4MuPass[i]);
	 	Hist->cut4ElPassHist->SetBinContent(i+1, cut4ElPass[i]);
	  	Hist->cut2L2LPassHist->SetBinContent(i+1, cut2L2LPass[i]);
		cut4MuPassW[i] = Hist->cut4MuPassHistW->GetBinContent(i+1);
		cut4ElPassW[i] = Hist->cut4ElPassHistW->GetBinContent(i+1);
		cut2L2LPassW[i] = Hist->cut2L2LPassHistW->GetBinContent(i+1);
	 }
}
void HiggsAnalysis::printMuonInfo()
{
	cout<<"-------------------"<<endl;
	cout<<"Muon Information"<<endl;
	cout<<"Staco"<<endl;	
	for(Int_t i = 0; i < event->mu_staco.n(); i++)
	{
		cout<<"Muon: "<< i<< " eta: "<< event->mu_staco[i].eta()<<" phi:" << event->mu_staco[i].phi();
		cout<<"  pT: "<< event->mu_staco[i].pt()<<endl;
		cout<<"Muon: "<< i<< " IDeta: "<< event->mu_staco[i].id_eta<<" IDphi:" << event->mu_staco[i].id_phi();
		cout<<"  IDpT: "<< event->mu_staco[i].id_pt<<endl;
		cout<<"Muon: "<< i<< " MSeta: "<< event->mu_staco[i].me_eta<<" MSphi:" << event->mu_staco[i].me_phi();
		cout<<"  MSpT: "<< event->mu_staco[i].me_pt<<endl;
	}
	cout<<"Calo"<<endl;	
	for(Int_t i = 0; i < event->mu_calo.n(); i++)
	{
		cout<<"Muon: "<< i<< " eta: "<< event->mu_calo[i].eta()<<" phi:" << event->mu_calo[i].phi();
		cout<<"  pT: "<< event->mu_calo[i].pt()<<endl;
		cout<<"Muon: "<< i<< " IDeta: "<< event->mu_calo[i].id_eta<<" IDphi:" << event->mu_calo[i].id_phi();
		cout<<"  IDpT: "<< event->mu_calo[i].id_pt<<endl;
		cout<<"Muon: "<< i<< " MSeta: "<< event->mu_calo[i].me_eta<<" MSphi:" << event->mu_calo[i].me_phi();
		cout<<"  MSpT: "<< event->mu_calo[i].me_pt<<endl;
	}
}
