#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <TString.h>
#include <TChain.h>
#include <TROOT.h>
#include "../MyAnalysis/macroDef.h"
#include <TStopwatch.h>
using namespace std;

// Reads tokens from a file delimited by '\n' and ','
string readFileName(TString filename) {
  vector<string> tokens;

  ifstream file(filename);
  if(!file.is_open()) return tokens;

  while(file.good()) {
    string line;
    getline(file, line);
    if(line == "") continue;

  	return line; 
  } // file.good()

}


// Reads tokens from a file delimited by '\n' and ','
vector<std::string> readAll(TString filename) {
  vector<string> tokens;

  ifstream file(filename);
  if(!file.is_open()) return tokens;

  while(file.good()) {
    string line;
    getline(file, line);
    if(line == "") continue;

    istringstream ss(line);
    while(ss.good()) {
      string token;
      getline(ss, token, ',');
      tokens.push_back(token);
    } // ss.good()
  } // file.good()

  return tokens;
}

// Returns all tokens that are root file names from input file
vector<string> getFilenames(TString filename) {
  vector<string> tokens = readAll(filename);
  vector<string> filenames;

  for(size_t i = 0; i < tokens.size(); i++) {
    if(tokens[i].find(".root") < tokens[i].length()) {
      cout << "Input file: " << tokens[i] << endl;
      filenames.push_back(tokens[i]);
    }
  }

  return filenames;
}

void runAnalysisGrid(Bool_t useNewGeo = false)
{
	
	cout<<"file: "<<readFileName("inputName.txt")<<endl;
	// Load the libraries that ROOT core compiled...
	gROOT->ProcessLine(".x $ROOTCOREDIR/scripts/load_packages.C");
	
	// Initializing the varibles....
	TString filePath;

	// For Keeping track of when to overwrite
	Bool_t fileOverWrite = true;

	// Actual analysis
	TString dataFileName ("input.txt");

	// To Count time	
	TStopwatch timet;
	timet.Start();
	TChain *phyData = new TChain ("physics");

	// Reading the file
	vector<string> inputFilenames = getFilenames(dataFileName);
	
	for(Int_t i = 0; i < inputFilenames.size(); i++)
	{
		phyData->Add(inputFilenames[i].c_str());
	}

	// Initialzing the cutFlow
	HiggsAnalysis *cutFlow = new HiggsAnalysis(phyData, true, readFileName("inputName.txt"), useNewGeo, doAnalysis::StdHZZllll);
	//HiggsAnalysis *cutFlow = new HiggsAnalysis(phyData, true, readFileName("inputName.txt"), useNewGeo, doAnalysis::Zllll);
	//HiggsAnalysis *cutFlow = new HiggsAnalysis(phyData, true, readFileName("inputName.txt"), useNewGeo, doAnalysis::trigeff4l);

	cutFlow->runningGrid = true;
	cutFlow->gridFileName = readFileName("inputName.txt");
	cutFlow->InitializeVar();
	// Setting up the printing
	cutFlow->SetupPrintEventList(fileOverWrite, "EventList");

	cutFlow->outputFilePath = "Output/EventSummary.root";
	

	// The call to the actual analysis
	int passedAllCut = 0;
	passedAllCut = cutFlow->AnalyzeTree();
	

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
		cout<<setw(16)<<cutFlow->cutElName[i]<<":\t"<<cutFlow->cutElPass[i]<<":\t"<<cutFlow->cutElLoosePass[i]<<endl;
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
