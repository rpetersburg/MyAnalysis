#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <TString.h>
#include <TChain.h>
#include <TROOT.h>
#include <TStopwatch.h>
#include "../MyAnalysis/macroDef.h"
using namespace std;
void runAnalysisRoot()
{
	// Load the libraries that ROOT core compiled...
	gROOT->ProcessLine(".x $ROOTCOREDIR/scripts/load_packages.C");
	
	Bool_t do2011_data 		= false;
	Bool_t do2011_ttH125 	= false;
	Bool_t do2011_VBF125 	= false;
	Bool_t do2011_WH125 	= false;
	Bool_t do2011_ZH125 	= false;
	
	Bool_t do2012_data 		= false;
	Bool_t do2012_ttH125 	= false;
	Bool_t do2012_VBF125 	= true;
	Bool_t do2012_ggH125 	= false;
	Bool_t do2012_WH125 	= false;
	Bool_t do2012_ZH125 	= false;
	Bool_t do2012_Spin	 	= false;
	
	Bool_t do2012_Zee 		= false;
	Bool_t do2012_Zmumu 	= false;

	// Initializing the varibles....
	TString filePath;

	// For Keeping track of when to overwrite
	Bool_t fileOverWrite = true;

	// Actual analysis
	ifstream dataFileName;
	if(do2011_data) 		dataFileName.open("DataFile/2011_data.txt");
	else if(do2011_ttH125) 	dataFileName.open("DataFile/2011_ttH125.txt");
	else if(do2011_VBF125) 	dataFileName.open("DataFile/2011_VBF125.txt");
	else if(do2011_WH125) 	dataFileName.open("DataFile/2011_WH125.txt");
	else if(do2011_ZH125) 	dataFileName.open("DataFile/2011_ZH125.txt");

	else if(do2012_data) 	dataFileName.open("DataFile/2012_data.txt");
	else if(do2012_ttH125) 	dataFileName.open("DataFile/2012_ttH125.txt");
	else if(do2012_VBF125) 	dataFileName.open("DataFile/2012_VBF125.txt");
	else if(do2012_ggH125) 	dataFileName.open("DataFile/2012_ggH125.txt");
	else if(do2012_WH125) 	dataFileName.open("DataFile/2012_WH125.txt");
	else if(do2012_ZH125) 	dataFileName.open("DataFile/2012_ZH125.txt");
	else if(do2012_Spin) 	dataFileName.open("DataFile/2012_Spin.txt");

	
	// To Count time	
	TStopwatch timet;
	timet.Start();
	TChain *phyData = new TChain ("physics");

	// Reading the TChain Files
	do{

 	  	dataFileName >> filePath;
 		// Checking if there is a file name or not
 		if (filePath.Length() <= 0 || filePath.Contains("#")) continue;
 		cout << filePath<<endl;
		// Adding it to the chain		
		phyData->Add(filePath.Data());
	}while (dataFileName.good());

	// Initialzing the cutFlow
	HiggsAnalysis *cutFlow = new HiggsAnalysis(phyData, false, "GridRunMismatch", true, doAnalysis::StdHZZllll);
	//HiggsAnalysis *cutFlow = new HiggsAnalysis(phyData, false, "GridRunMismatch", true, doAnalysis::trigeff4l);
	//HiggsAnalysis *cutFlow = new HiggsAnalysis(phyData, false, "GridRunMismatch", true, doAnalysis::Zllll);
	cutFlow->InitializeVar();
	// Setting up the printing
	cutFlow->SetupPrintEventList(fileOverWrite, "EventList");
	if(do2011_data) 				cutFlow->outputFilePath = "Output/mc11d_data.root";
	else if(do2011_ttH125) 			cutFlow->outputFilePath = "Output/mc11d_ttH125.root";
	else if(do2011_VBF125) 			cutFlow->outputFilePath = "Output/mc11d_VBFH125.root";
	else if(do2011_WH125) 			cutFlow->outputFilePath = "Output/mc11d_WH125.root";
	else if(do2011_ZH125) 			cutFlow->outputFilePath = "Output/mc11d_ZH125.root";

	else if(do2012_data) 			cutFlow->outputFilePath = "Output/mc12c_data.root";
	else if(do2012_ttH125) 			cutFlow->outputFilePath = "Output/mc12c_ttH125.root";
	else if(do2012_VBF125) 			cutFlow->outputFilePath = "Output/mc12c_VBFH125.root";
	else if(do2012_ggH125) 			cutFlow->outputFilePath = "Output/mc12c_ggH125.root";
	else if(do2012_WH125) 			cutFlow->outputFilePath = "Output/mc12c_WH125.root";
	else if(do2012_ZH125) 			cutFlow->outputFilePath = "Output/mc12c_ZH125.root";
	else if(do2012_Spin) 			cutFlow->outputFilePath = "Output/mc12c_Spin.root";

	// The call to the actual analysis
	int passedAllCut = 0;
	passedAllCut = cutFlow->AnalyzeTree();
	
	// MC Studies
	//cutFlow->printMCInfo(55576);

	// For Debug
	
	//cutFlow->printDebug(37543, false);
	//cutFlow->printDebug(26077, false);
	//cutFlow->printDebug(34624, false);
	//cutFlow->printDebug(30205, false);
	//cutFlow->printDebug(30247, true);
	//cutFlow->printDebug(5604);
	//cutFlow->printDebug(8185);
	//cutFlow->printDebug(42353);
	//cutFlow->printDebug(1968);
	//cutFlow->printDebug(24605);
	//cutFlow->printMCInfo(58);

	// Saving the final Histrogram
	cutFlow->SaveHist(fileOverWrite);	
	fileOverWrite = false;
	// Stop the time
	timet.Stop();

	// Final Output
	cout<<"Passed Cut: "<< passedAllCut<<endl;
	cout<<"Time (Real): "<<timet.RealTime()<<" Time (CPU):"<<timet.CpuTime()<<endl;
	cout<<"Total Cut Flow"<<endl; 
	
	// Printing the Final Number
	for(Int_t i = 0; i < cutFlow->nCut; i++)
	{
		cout<<setw(16)<<cutFlow->cutName[i]<<":\t"<<cutFlow->cutPass[i]<<endl;
	}
	// Printing the muon final Number
	cout <<endl<<"Muon Cut Flow"<<endl; 
	for(Int_t i = 0; i < cutFlow->nMuCut; i++)
	{
		cout<<setw(16)<<cutFlow->cutMuName[i]<<":\t"<<cutFlow->cutMuPass[i]<<endl;
	}
	// Printing the electron final Number
	cout <<endl<<"Electron Cut Flow"<<endl; 
	for(Int_t i = 0; i < cutFlow->nElCut; i++)
	{
		cout<<setw(16)<<cutFlow->cutElName[i]<<":\t"<<cutFlow->cutElPass[i]<<"\t:\t"<<cutFlow->cutElLoosePass[i]<<endl;
	}
	// Printing the jets final Number
	cout <<endl<<"Jets Cut Flow"<<endl; 
	for(Int_t i = 0; i < cutFlow->nJetsCut; i++)
	{
		cout<<setw(16)<<cutFlow->cutJetsName[i]<<":\t"<<cutFlow->cutJetsPass[i]<<endl;
	}
	// Printing the channel specific Number
	cout <<endl<<"Channel Flow\t\t4MU\t\t4E\t\t2L2L"<<endl; 
	for(Int_t i = 0; i < cutFlow->nCH; i++)
	{
		cout<<setw(12)<<cutFlow->cutCHName[i]<<":\t\t";
		cout<<setw(7)<<cutFlow->cut4MuPass[i]<<"\t\t";
		cout<<setw(7)<<cutFlow->cut4ElPass[i]<<"\t\t";
		cout<<setw(7)<<cutFlow->cut2L2LPass[i]<<"\t\t";
		cout<<setw(7)<<cutFlow->cutlleePass[i]<<endl;
		
	}
	if(cutFlow->doWeight)
	{
		cout<<endl<<"Cut Flow with Weights"<<endl;

		// Printing the Final Number
		for(Int_t i = 0; i < cutFlow->nCut; i++)
		{
			cout<<setw(16)<<setprecision(7)<<cutFlow->cutName[i]<<":\t"<<cutFlow->cutPassW[i]<<endl;
		}
		// Printing the channel specific Number
		cout <<endl<<"Channel Flow\t\t4MU\t\t4E\t\t2L2L"<<endl; 
		for(Int_t i = 0; i < cutFlow->nCH; i++)
		{
			cout<<setw(12)<<cutFlow->cutCHName[i]<<":\t\t";
			cout<<setw(7)<<setprecision(6)<<cutFlow->cut4MuPassW[i]<<"\t\t";
			cout<<setw(7)<<setprecision(6)<<cutFlow->cut4ElPassW[i]<<"\t\t";
			cout<<setw(7)<<setprecision(6)<<cutFlow->cut2L2LPassW[i]<<endl;
			
		}
	}

	// Printing the production Channel Number
	cout <<endl<<"Production Channel\t\t4MU\t\t4E\t\t2L2L"<<endl; 
	for(Int_t i = 0; i < cutFlow->nProdCH; i++)
	{
		cout<<setw(12)<<cutFlow->prodCHName[i]<<":\t\t";
		cout<<setw(7)<<cutFlow->prodCH4Mu[i]<<"\t\t";
		cout<<setw(7)<<cutFlow->prodCH4El[i]<<"\t\t";
		cout<<setw(7)<<cutFlow->prodCH2L2L[i]<<"\t\t";
		cout<<setw(7)<<cutFlow->prodCHllee[i]<<endl;
		
	}

	cout <<endl<<"Truth quad type"<<endl; 
	for(Int_t i = 0; i < cutFlow->nListTruthQuadType; i++)
	{
		cout<<setw(8)<<cutFlow->truthQuadType[i]<<":\t\t";
		cout<<setw(5)<<cutFlow->nTruthQuadType[i]<<endl;
		
	}
	// Deleting the vars
	delete cutFlow;
	delete phyData;


}
