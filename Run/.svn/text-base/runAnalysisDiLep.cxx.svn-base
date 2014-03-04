#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <TString.h>
#include <TChain.h>
#include <TROOT.h>
#include <TStopwatch.h>

using namespace std;
void runAnalysisDiLep()
{
	// Load the libraries that ROOT core compiled...
	gROOT->ProcessLine(".x $ROOTCOREDIR/scripts/load_packages.C");
	
	Bool_t doZee 	= true;
	Bool_t doZmumu	= false;
	// Initializing the varibles....
	TString filePath;

	// For Keeping track of when to overwrite
	Bool_t fileOverWrite = true;

	// Actual analysis
	ifstream dataFileName;
	if(doZee)  			dataFileName.open("DataFile/DataFileZee.txt");	
	else if(doZmumu)dataFileName.open("DataFile/DataFileZmumu.txt");	
		
	// To Count time	
	TStopwatch timet;
	timet.Start();
	TChain *phyData = new TChain ("physics");

	// Reading the TChain Files
	do{

 	  	dataFileName >> filePath;
 		// Checking if there is a file name or not
 		if (filePath.Length() <= 0) continue;
 		cout << filePath<<endl;
		// Adding it to the chain		
		phyData->Add(filePath.Data());
	}while (dataFileName.good());

	// Initialzing the cutFlow
	DiLepAnalysis *cutFlow = new DiLepAnalysis(phyData);
	cutFlow->InitializeVar();
	// Setting up the printing
	//if(doMC2011 || doMC2012) cutFlow->SetupPrintEventList(fileOverWrite, "EventList");
	//else if(doData2011 || doData2012) cutFlow->SetupPrintEventList(fileOverWrite, "EventListData");
	//else if(doggF2011 || doBkg2012 || doMCTest) cutFlow->SetupPrintEventList(fileOverWrite, "EventListTest");
	//else  cutFlow->SetupPrintEventList(fileOverWrite, "EventListTest1");
	cutFlow->SetupPrintEventList(fileOverWrite, "EventList");
	
	if(doZee) 			cutFlow->outputFilePath = "Output/mc12_Zee.root";
	else if(doZmumu)	cutFlow->outputFilePath = "Output/mc12_Zmumu.root";
	
	// The call to the actual analysis
	int passedAllCut = 0;
	passedAllCut = cutFlow->AnalyzeTree();
	
	// MC Studies
	//cutFlow->printMCInfo(30028);

	// For Debug
	//cutFlow->printDebug(92658);
	//cutFlow->printDebug(2669);
	//cutFlow->printDebug(30128);
	//cutFlow->printDebug(30153);
	//cutFlow->printDebug(30609);
	//cutFlow->printDebug(1968);
	//cutFlow->printDebug(24605);
	//cutFlow->printMCInfo(24605);

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
		cout<<setw(16)<<cutFlow->cutElName[i]<<":\t"<<cutFlow->cutElPass[i]<<endl;
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
		cout<<setw(7)<<cutFlow->cut2L2LPass[i]<<endl;
		
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
		cout<<setw(7)<<cutFlow->prodCH2L2L[i]<<endl;
		
	}


	// Deleting the vars
	delete cutFlow;
	delete phyData;


}
