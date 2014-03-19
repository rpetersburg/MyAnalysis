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
	
	Bool_t doMC2011 		= false;
	Bool_t doMC2011NewGeo	= false;
	Bool_t doMC2012 		= false;
	Bool_t doMC2012NewGeo	= true;
	Bool_t doData2011 		= false;
	Bool_t doData2012 		= false;
	Bool_t doggF2011 		= false;
	Bool_t doBkg2012 		= false;
	Bool_t doMCTest 		= false;
	Bool_t doMix1 			= false;
	Bool_t doMix2 			= false;
	Bool_t doMix3 			= false;
	Bool_t doSpin 			= false;
	Bool_t doSpin2011 		= false;
	Bool_t doWH2011 		= false;
	Bool_t doWH2012 		= false;
	Bool_t doWH2012NewGeo	= false;
	Bool_t doZH2012 		= false;
	Bool_t doZH2012NewGeo	= false;
	Bool_t doMCatNlo 		= false;
	Bool_t dottH2011 		= false;
	Bool_t dottH2012 		= false;
	Bool_t dottH2012NewGeo	= false;
	Bool_t doggH121 		= false;
	Bool_t doggH123p5 		= false;
	Bool_t dodataTest2012 	= false;
	Bool_t dodataTest2011 	= false;

	// Initializing the variables....
	TString filePath;

	// For Keeping track of when to overwrite
	Bool_t fileOverWrite = true;

	// Actual analysis
	ifstream dataFileName;
	if(doMC2011) 				dataFileName.open("DataFile/DataFileMC2011Chain.txt");
	else if(doMC2011NewGeo)		dataFileName.open("DataFile/DataFileMC2011_NewGeo_Chain.txt");
	else if(doMC2012) 			dataFileName.open("DataFile/DataFileMC2012Chain.txt");
	else if(doMC2012NewGeo)  	dataFileName.open("DataFile/DataFileMC2012_NewGeo_Chain.txt");
	else if(doData2011) 		dataFileName.open("DataFile/DataFileData2011Chain.txt");
	else if(doData2012)  		dataFileName.open("DataFile/DataFileData2012Chain.txt");
	else if(doggF2011)  		dataFileName.open("DataFile/DataFileData2011ggF.txt");	
	else if(doBkg2012)  		dataFileName.open("DataFile/DataFile2012Bkg.txt");	
	else if(doMCTest)  			dataFileName.open("DataFile/DataFileMCtest.txt");	
	else if(doMix1)  			dataFileName.open("DataFile/DataFile2012Mix1.txt");	
	else if(doMix2) 			dataFileName.open("DataFile/DataFile2012Mix2.txt");	
	else if(doMix3)  			dataFileName.open("DataFile/DataFile2012Mix3.txt");	
	else if(doSpin2011)  		dataFileName.open("DataFile/DataFile2011Spin.txt");	
	else if(doSpin)  			dataFileName.open("DataFile/DataFile2012Spin.txt");	
	else if(doWH2011)  			dataFileName.open("DataFile/DataFileWHMC2011Chain.txt");	
	else if(doWH2012)  			dataFileName.open("DataFile/DataFileWHMC2012Chain.txt");	
	else if(doWH2012NewGeo) 	dataFileName.open("DataFile/DataFileWHMC2012_NewGeo_Chain.txt");	
	else if(doZH2012)  			dataFileName.open("DataFile/DataFileZHMC2012Chain.txt");	
	else if(doZH2012NewGeo)  	dataFileName.open("DataFile/DataFileZHMC2012_NewGeo_Chain.txt");	
	else if(doMCatNlo)  		dataFileName.open("DataFile/DataFileMCAtNLO2012Chain.txt");	
	else if(dottH2011)  		dataFileName.open("DataFile/DataFilettHMC2011Chain.txt");	
	else if(dottH2012)  		dataFileName.open("DataFile/DataFilettHMC2012Chain.txt");	
	else if(dottH2012NewGeo)  	dataFileName.open("DataFile/DataFilettHMC2012_NewGeo_Chain.txt");	
	else if(doggH121)  			dataFileName.open("DataFile/DataFileggH121MC2012Chain.txt");	
	else if(doggH123p5)  		dataFileName.open("DataFile/DataFileggH123p5MC2012Chain.txt");	
	else if(dodataTest2012)  	dataFileName.open("DataFile/DataFileDataTest2012.txt");	
	else if(dodataTest2011)  	dataFileName.open("DataFile/DataFileDataTest2011.txt");	
	
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
	HiggsAnalysis *cutFlow = new HiggsAnalysis(phyData, false, "GridRunMismatch", true, doAnalysis::StdHZZllll);
	//HiggsAnalysis *cutFlow = new HiggsAnalysis(phyData, false, "GridRunMismatch", true, doAnalysis::trigeff4l);
	//HiggsAnalysis *cutFlow = new HiggsAnalysis(phyData, false, "GridRunMismatch", true, doAnalysis::Zllll);
	cutFlow->InitializeVar();
	// Setting up the printing
	cutFlow->SetupPrintEventList(fileOverWrite, "EventList");
	
	if(doMC2011) 				cutFlow->outputFilePath = "Output/mc11a_VBFH125.root";
	else if(doMC2011NewGeo)		cutFlow->outputFilePath = "Output/mc11d_VBFH125.root";
	else if(doMC2012) 			cutFlow->outputFilePath = "Output/mc12a_VBFH125.root";
	else if(doMC2012NewGeo) 	cutFlow->outputFilePath = "Output/mc12c_VBFH125.root";
	else if(doData2011) 		cutFlow->outputFilePath = "Output/data11.root";
	else if(doData2012) 		cutFlow->outputFilePath = "Output/data12.root";
	else if(doggF2011) 			cutFlow->outputFilePath = "Output/mc11a_ggH126.root";
	else if(doBkg2012) 			cutFlow->outputFilePath = "Output/mc12a_Bkg.root";
	else if(doMCTest) 			cutFlow->outputFilePath = "Output/mcTest.root";
	else if(doMix1) 			cutFlow->outputFilePath = "Output/mc12a_Mix_1.root";
	else if(doMix2) 			cutFlow->outputFilePath = "Output/mc12a_Mix_2.root";
	else if(doMix3) 			cutFlow->outputFilePath = "Output/mc12a_Mix_3.root";
	else if(doSpin2011) 		cutFlow->outputFilePath = "Output/mc11a_spinJHU.root";
	else if(doSpin) 			cutFlow->outputFilePath = "Output/mc12a_spinJHU.root";
	else if(doWH2011) 			cutFlow->outputFilePath = "Output/mc11a_WH125.root";
	else if(doWH2012) 			cutFlow->outputFilePath = "Output/mc12a_WH125.root";
	else if(doWH2012NewGeo) 	cutFlow->outputFilePath = "Output/mc12c_WH125.root";
	else if(doZH2012) 			cutFlow->outputFilePath = "Output/mc12a_ZH125.root";
	else if(doZH2012NewGeo) 	cutFlow->outputFilePath = "Output/mc12c_ZH125.root";
	else if(doMCatNlo) 			cutFlow->outputFilePath = "Output/mc12a_MCatNLO.root";
	else if(dottH2011) 			cutFlow->outputFilePath = "Output/mc11a_ttH125.root";
	else if(dottH2012) 			cutFlow->outputFilePath = "Output/mc12a_ttH125.root";
	else if(dottH2012NewGeo)	cutFlow->outputFilePath = "Output/mc12c_ttH125.root";
	else if(doggH121) 			cutFlow->outputFilePath = "Output/mc12a_ggH121.root";
	else if(doggH123p5) 		cutFlow->outputFilePath = "Output/mc12a_ggH123p5.root";
	else if(dodataTest2012) 	cutFlow->outputFilePath = "Output/dataTest12.root";
	else if(dodataTest2011) 	cutFlow->outputFilePath = "Output/dataTest11.root";

	// The call to the actual analysis
	int passedAllCut = 0;
	passedAllCut = cutFlow->AnalyzeTree();
	
	// MC Studies
	//cutFlow->printMCInfo(30028);

	// For Debug
	
	// cutFlow->printDebug(45623, false);
	// cutFlow->printDebug(29461, false);
	// cutFlow->printDebug(393, false);
	// cutFlow->printDebug(32632, false);
	// cutFlow->printDebug(24565, true);
	// cutFlow->printDebug(5604);
	// cutFlow->printDebug(8185);
	// cutFlow->printDebug(42353);
	// cutFlow->printDebug(1968);
	// cutFlow->printDebug(24605);
	// cutFlow->printMCInfo(42353);

	// // Saving the final Histogram
	// cutFlow->SaveHist(fileOverWrite);	
	// fileOverWrite = false;
	// // Stop the time
	// timet.Stop();

	// // Final Output
	// cout<<"Passed Cut: "<< passedAllCut<<endl;
	// cout<<"Time (Real): "<<timet.RealTime()<<" Time (CPU):"<<timet.CpuTime()<<endl;
	// cout<<"Total Cut Flow"<<endl; 
	
	// // Printing the Final Number
	// for(Int_t i = 0; i < cutFlow->nCut; i++)
	// {
		// cout<<setw(16)<<cutFlow->cutName[i]<<":\t"<<cutFlow->cutPass[i]<<endl;
	// }
	// // Printing the muon final Number
	// cout <<endl<<"Muon Cut Flow"<<endl; 
	// for(Int_t i = 0; i < cutFlow->nMuCut; i++)
	// {
		// cout<<setw(16)<<cutFlow->cutMuName[i]<<":\t"<<cutFlow->cutMuPass[i]<<endl;
	// }
	// // Printing the electron final Number
	// cout <<endl<<"Electron Cut Flow"<<endl; 
	// for(Int_t i = 0; i < cutFlow->nElCut; i++)
	// {
		// cout<<setw(16)<<cutFlow->cutElName[i]<<":\t"<<cutFlow->cutElPass[i]<<"\t:\t"<<cutFlow->cutElLoosePass[i]<<endl;
	// }
	// // Printing the jets final Number
	// cout <<endl<<"Jets Cut Flow"<<endl; 
	// for(Int_t i = 0; i < cutFlow->nJetsCut; i++)
	// {
		// cout<<setw(16)<<cutFlow->cutJetsName[i]<<":\t"<<cutFlow->cutJetsPass[i]<<endl;
	// }
	// // Printing the channel specific Number
	// cout <<endl<<"Channel Flow\t\t4MU\t\t4E\t\t2L2L"<<endl; 
	// for(Int_t i = 0; i < cutFlow->nCH; i++)
	// {
		// cout<<setw(12)<<cutFlow->cutCHName[i]<<":\t\t";
		// cout<<setw(7)<<cutFlow->cut4MuPass[i]<<"\t\t";
		// cout<<setw(7)<<cutFlow->cut4ElPass[i]<<"\t\t";
		// cout<<setw(7)<<cutFlow->cut2L2LPass[i]<<"\t\t";
		// cout<<setw(7)<<cutFlow->cutlleePass[i]<<endl;
		
	// }
	// if(cutFlow->doWeight)
	// {
		// cout<<endl<<"Cut Flow with Weights"<<endl;

		// // Printing the Final Number
		// for(Int_t i = 0; i < cutFlow->nCut; i++)
		// {
			// cout<<setw(16)<<setprecision(7)<<cutFlow->cutName[i]<<":\t"<<cutFlow->cutPassW[i]<<endl;
		// }
		// // Printing the channel specific Number
		// cout <<endl<<"Channel Flow\t\t4MU\t\t4E\t\t2L2L"<<endl; 
		// for(Int_t i = 0; i < cutFlow->nCH; i++)
		// {
			// cout<<setw(12)<<cutFlow->cutCHName[i]<<":\t\t";
			// cout<<setw(7)<<setprecision(6)<<cutFlow->cut4MuPassW[i]<<"\t\t";
			// cout<<setw(7)<<setprecision(6)<<cutFlow->cut4ElPassW[i]<<"\t\t";
			// cout<<setw(7)<<setprecision(6)<<cutFlow->cut2L2LPassW[i]<<endl;
			
		// }
	// }

	// // Printing the production Channel Number
	// cout <<endl<<"Production Channel\t\t4MU\t\t4E\t\t2L2L"<<endl; 
	// for(Int_t i = 0; i < cutFlow->nProdCH; i++)
	// {
		// cout<<setw(12)<<cutFlow->prodCHName[i]<<":\t\t";
		// cout<<setw(7)<<cutFlow->prodCH4Mu[i]<<"\t\t";
		// cout<<setw(7)<<cutFlow->prodCH4El[i]<<"\t\t";
		// cout<<setw(7)<<cutFlow->prodCH2L2L[i]<<"\t\t";
		// cout<<setw(7)<<cutFlow->prodCHllee[i]<<endl;
		
	// }

	// cout <<endl<<"Truth quad type"<<endl; 
	// for(Int_t i = 0; i < cutFlow->nListTruthQuadType; i++)
	// {
		// cout<<setw(8)<<cutFlow->truthQuadType[i]<<":\t\t";
		// cout<<setw(5)<<cutFlow->nTruthQuadType[i]<<endl;
		
	// }
	// Deleting the vars
	delete cutFlow;
	delete phyData;


}
