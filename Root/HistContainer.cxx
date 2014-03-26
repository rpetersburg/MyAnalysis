#include <stdlib.h>
#include <string>
#include "MyAnalysis/HistContainer.h"


using namespace std;
////////////////////////////////////////////////////////////////////////////////////////
//							Constructor and Destructor
////////////////////////////////////////////////////////////////////////////////////////

HistContainer::HistContainer(Int_t nCut, Int_t nMuCut, Int_t nElCut, Int_t nJetsCut, Int_t nCH)
{
	// Weight 
	weight = 1;

	// Histograms. Need to create a functions for them	
	hist4MuMUnconstrained = new TH1D("hist4MuMUnconstrained", "4Mu QuadLepton Mass;Energy [MeV]; Events", 408, 80*1000, 1000*1000);
	hist4ElMUnconstrained = new TH1D("hist4ElMUnconstrained", "4El QuadLepton Mass;Energy [MeV]; Events", 408, 80*1000, 1000*1000);
	hist2El2MuMUnconstrained = new TH1D("hist2El2MuMUnconstrained", "2El2Mu QuadLepton Mass;Energy [MeV]; Events", 408, 80*1000, 1000*1000);
	hist2Mu2ElMUnconstrained = new TH1D("hist2Mu2ElMUnconstrained", "2Mu2El QuadLepton Mass;Energy [MeV]; Events", 408, 80*1000, 1000*1000);
	
	hist4MuMFSR = new TH1D("hist4MuMFSR", "4Mu QuadLepton Mass;Energy [MeV]; Events", 408, 80*1000, 1000*1000);
	hist4ElMFSR = new TH1D("hist4ElMFSR", "4El QuadLepton Mass;Energy [MeV]; Events", 408, 80*1000, 1000*1000);
	hist2El2MuMFSR = new TH1D("hist2El2MuMFSR", "2El2Mu QuadLepton Mass;Energy [MeV]; Events", 408, 80*1000, 1000*1000);
	hist2Mu2ElMFSR = new TH1D("hist2Mu2ElMFSR", "2Mu2El QuadLepton Mass;Energy [MeV]; Events", 408, 80*1000, 1000*1000);
	
	hist4MuMConstrained = new TH1D("hist4MuMConstrained", "4Mu QuadLepton Mass;Energy [MeV]; Events", 408, 80*1000, 1000*1000);
	hist4ElMConstrained = new TH1D("hist4ElMConstrained", "4El QuadLepton Mass;Energy [MeV]; Events", 408, 80*1000, 1000*1000);
	hist2El2MuMConstrained = new TH1D("hist2El2MuMConstrained", "2El2Mu QuadLepton Mass;Energy [MeV]; Events", 408, 80*1000, 1000*1000);
	hist2Mu2ElMConstrained = new TH1D("hist2Mu2ElMConstrained", "2Mu2El QuadLepton Mass;Energy [MeV]; Events", 408, 80*1000, 1000*1000);
		
	// Counting Histrograms
	cutPassHist = new TH1D("cutPassHist", "cutPassHist;CutNum;Events", nCut, 0, nCut);
	cutMuPassHist = new TH1D("cutMuPassHist", "cutMuPassHist;CutNum;Events", nMuCut, 0, nMuCut);
	cutElPassHist = new TH1D("cutElPassHist", "cutElPassHist;CutNum;Events", nElCut, 0, nElCut);
	cutJetsPassHist = new TH1D("cutJetsPassHist", "cutJetsPassHist;CutNum;Events", nJetsCut, 0, nJetsCut);
	cut4MuPassHist = new TH1D("cut4MuPassHist", "cut4MuPassHist;CutNum;Events", nCH, 0, nCH);
	cut4ElPassHist = new TH1D("cut4ElPassHist", "cut4ElPassHist;CutNum;Events", nCH, 0, nCH);
	cut2L2LPassHist = new TH1D("cut2L2LPassHist", "cut2L2LPassHist;CutNum;Events", nCH, 0, nCH);
	truthQuadHist = new TH1D("truthQuadHist", "truthQuadHist;type;Events", 10, 0, 10);
	// Weights
	cutPassHistW = new TH1D("cutPassHistWeight", "cutPassHistWeight;CutNum;Events", nCut, 0, nCut);
	cut4MuPassHistW = new TH1D("cut4MuPassHistWeight", "cut4MuPassHistWeight;CutNum;Events", nCH, 0, nCH);
	cut4ElPassHistW = new TH1D("cut4ElPassHistWeight", "cut4ElPassHistWeight;CutNum;Events", nCH, 0, nCH);
	cut2L2LPassHistW = new TH1D("cut2L2LPassHistWeight", "cut2L2LPassHistWeight;CutNum;Events", nCH, 0, nCH);
	cutlleePassHistW = new TH1D("cutlleePassHistWeight", "cutlleePassHistWeight;CutNum;Events", nCH, 0, nCH);

	// Muon Histrogram
	muAuthorHist = new TH1D *[3];
	muPTHist = new TH1D *[3];
	muEtaHist = new TH1D *[3];
	muD0Hist = new TH1D *[3];
	muZ0Hist = new TH1D *[3];
	muOverlapHist = new TH1D *[3];
	TString muName [3] = {"Staco", "Calo", "StandAlone"};
	for(Int_t i = 0; i < 3; i++)
	{
		muAuthorHist[i] = new TH1D("authorHist"+muName[i], "authorHist"+muName[i]+";Author;Event", 20, 0, 20);
		muPTHist[i] = new TH1D("pTHist"+muName[i], "pTHist"+muName[i]+";pT[Mev];Event", 100, 0, 250*1000);
		muEtaHist[i] = new TH1D("etaHist"+muName[i], "etaHist"+muName[i]+";eta;Event", 50, -4, 4);
		muD0Hist[i] = new TH1D("d0Hist"+muName[i], "d0Hist"+muName[i]+";d0[mm];Event", 100, -5, 5);
		muZ0Hist[i] = new TH1D("z0Hist"+muName[i], "z0Hist"+muName[i]+";z0[mm];Event", 100, -20, 20);
		muOverlapHist[i] = new TH1D("overlapHist"+muName[i], "overlapHist"+muName[i]+";DeltaR;Event", 50, 0, 5);
	}

	// Electron Histrogram
	elAuthorHist = new TH1D("authorHistElGFS", "authorHistElGFS;Author;Event", 20, 0, 20);
	elETHist = new TH1D("ETHistElGFS", "ETHistElGFS;pT[Mev];Event", 100, 0, 250*1000);;
	elEtaHist = new TH1D("etaHistElGFS", "etaHistElGFS;eta;Event", 50, -4, 4);
	elZ0Hist = new TH1D("z0HistElGFS", "z0HistElGFS;z0[mm];Event", 100, -20, 20);
	
	// QuadLepton
	quadInitPtHist = new TH1D **[4];
	quadFinalPtHist = new TH1D **[4];
	Z1MassHist = new TH1D *[4];
	Z2M4lHist = new TH2D *[4];
	JPsiHist = new TH1D *[4];
	TrackIsoHist = new TH1D *[4];
	muCaloIsoHist = new TH1D *[4];
	elCaloIsoHist = new TH1D *[4];
	muSigHist = new TH1D *[4];
	elSigHist = new TH1D *[4];
	TString chName [4] = {"4Mu", "4El", "2El2Mu", "2Mu2El"};
	TString pTName [4] = {"1", "2", "3", "4"};
	
	for(Int_t i = 0; i < 4; i++)
	{
		quadInitPtHist[i] = new TH1D *[4];
		quadFinalPtHist[i] = new TH1D *[4];
		for(Int_t j = 0; j < 4; j++)
		{
			quadInitPtHist[i][j] = new TH1D("pT"+pTName[j]+"InitHist"+chName[i], "pT"+pTName[j]+"InitHist"+chName[i]+";pT[Mev];Event", 
					100, 0, 250*1000);
			quadFinalPtHist[i][j] = new TH1D("pT"+pTName[j]+"FinalHist"+chName[i], "pT"+pTName[j]+"FinalHist"+chName[i]+"pT[Mev];Event", 
					100, 0, 250*1000);
		}
		
		Z1MassHist[i] = new TH1D("z1MassHist"+chName[i], "z1MassHist"+chName[i]+"MZ1[Mev];Event", 100, 0, 200*1000) ;
		Z2M4lHist[i] = new TH2D("z2MassHist"+chName[i], "z1MassHist"+chName[i]+";MZ2[Mev];M4l[Mev]", 
				100, 0, 100*1000, 100, 0 ,250*1000 );
		JPsiHist[i] = new TH1D("JPsiHist"+chName[i], "JPsiHist"+chName[i]+";MCrossPair[Mev];Event", 100, 0, 200*1000);
		TrackIsoHist[i] = new TH1D("TrackIsoHist"+chName[i], "TrackIsoHist"+chName[i]+";TrackIso;Event", 100, -2, 2);
		muCaloIsoHist[i] = new TH1D("muCaloIsoHist"+chName[i], "muCaloIsoHist"+chName[i]+";muCaloIsoHist;Event", 100, -2, 2);
		elCaloIsoHist[i] = new TH1D("elCaloIsoHist"+chName[i], "elCaloIsoHist"+chName[i]+";elCaloIsoHist;Event", 100, -2, 2);
		muSigHist[i] = new TH1D("muSigHist"+chName[i], "muSigHist"+chName[i]+";muImpactSig;Event", 100, -2, 10);
		elSigHist[i] = new TH1D("elSigHist"+chName[i], "elSigHist"+chName[i]+";elImpactSig;Event", 100, -2, 10);

	}

	muInteractionHist = new TH1D("muInteraction", "AverageCrossing;mu;Event", 50, 0, 50);
	muInteractionHistPileupW = new TH1D("muInteractionPileupW", "AverageCrossing;mu;Event", 50, 0, 50);

}

HistContainer::~HistContainer()
{
		// Gives a memory error
	//delete cutName;
	delete hist4MuMUnconstrained;
	delete hist4ElMUnconstrained;
	delete hist2El2MuMUnconstrained;
	delete hist2Mu2ElMUnconstrained;	

	delete hist4MuMFSR;
	delete hist4ElMFSR;
	delete hist2El2MuMFSR;
	delete hist2Mu2ElMFSR;
	
	delete hist4MuMConstrained;
	delete hist4ElMConstrained;
	delete hist2El2MuMConstrained;
	delete hist2Mu2ElMConstrained;
	
	delete cutPassHist;
	delete cutMuPassHist;
	delete cutElPassHist;
	delete cutJetsPassHist;
	delete cut4MuPassHist;
	delete cut4ElPassHist;
	delete cut2L2LPassHist;
	
	delete cutPassHistW;
	delete cut4MuPassHistW;
	delete cut4ElPassHistW;
	delete cut2L2LPassHistW;
	
	delete elAuthorHist;
	delete elETHist;
	delete elEtaHist;
	delete elZ0Hist;

	for(Int_t i = 0; i < 3; i++)
	{
		delete muAuthorHist[i];
		delete muPTHist[i];
		delete muEtaHist[i];
		delete muD0Hist[i];
		delete muZ0Hist[i];
		delete muOverlapHist[i];
	}
	delete[] muAuthorHist;
	delete[] muPTHist;
	delete[] muEtaHist;
	delete[] muD0Hist;
	delete[] muZ0Hist;
	delete[] muOverlapHist;

	for(Int_t i = 0; i < 4; i++)
	{
		for(Int_t j = 0; j < 4; j++)
		{
			delete	quadInitPtHist[i][j];
			delete quadFinalPtHist[i][j];
		}

		delete[] quadInitPtHist[i];
		delete[] quadFinalPtHist[i];
		
		delete Z1MassHist[i];
		delete Z2M4lHist[i];
		delete JPsiHist[i];
		delete TrackIsoHist[i];
		delete muCaloIsoHist[i];
		delete elCaloIsoHist[i];
		delete muSigHist[i];
		delete elSigHist[i];
	}
	delete[] quadInitPtHist;
	delete[] quadFinalPtHist;
	delete[] Z1MassHist;
	delete[] Z2M4lHist;
	delete[] JPsiHist;
	delete[] TrackIsoHist;
	delete[] muCaloIsoHist;
	delete[] elCaloIsoHist;
	delete[] muSigHist;
	delete[] elSigHist;
	delete muInteractionHist;
	delete muInteractionHistPileupW;

	delete truthQuadHist;
}

void HistContainer::SaveHist(TFile *output, TString fileNamePart, TH1F* countingHist)
{
	
	truthQuadHist->Write();
	// Saving the Counting Histrograms
	output->mkdir("Counting");
	output->cd("Counting");
	cutPassHist->Write();
	cutMuPassHist->Write();
	cutElPassHist->Write();
	cutJetsPassHist->Write();
	cut4MuPassHist->Write();
	cut4ElPassHist->Write();
	cut2L2LPassHist->Write();
	cutPassHistW->Write();
	cut4MuPassHistW->Write();
	cut4ElPassHistW->Write();
	cut2L2LPassHistW->Write();
	cutlleePassHistW->Write();

	// Muon Hist
	output->cd();
	output->mkdir("Muon");
	output->cd("Muon");
	muAuthorHist[0]->Write();
	muAuthorHist[1]->Write();
	muAuthorHist[2]->Write();
	muPTHist[0]->Write();
	muPTHist[1]->Write();
	muPTHist[2]->Write();
	muEtaHist[0]->Write();
	muEtaHist[1]->Write();
	muEtaHist[2]->Write();
	muD0Hist[0]->Write();
	muD0Hist[1]->Write();
	muZ0Hist[0]->Write();
	muZ0Hist[1]->Write();
	muOverlapHist[1]->Write();
	muOverlapHist[2]->Write();

	// Electron hist
	output->cd(); 
	output->mkdir("Electron");
	output->cd("Electron");
	elAuthorHist->Write();
	elETHist->Write();
	elEtaHist->Write();
	elZ0Hist->Write();

	// quadLepton hist
	TString chName[4] = {"4Mu", "4El", "2El2Mu", "2Mu2El"};
	for(Int_t i = 0; i < 4; i++)
	{
		output->cd();		
		output->mkdir(chName[i]);
		output->cd(chName[i]);
		quadInitPtHist[i][0]->Write();	
		quadInitPtHist[i][1]->Write();	
		quadInitPtHist[i][2]->Write();	
		quadInitPtHist[i][3]->Write();
		quadFinalPtHist[i][0]->Write();	
		quadFinalPtHist[i][1]->Write();	
		quadFinalPtHist[i][2]->Write();	
		quadFinalPtHist[i][3]->Write();

		Z1MassHist[i]->Write();
		Z2M4lHist[i]->Write();
	 	JPsiHist[i]->Write();
	 	TrackIsoHist[i]->Write();
	 	muCaloIsoHist[i]->Write();
	 	elCaloIsoHist[i]->Write();
	 	muSigHist[i]->Write();
	 	elSigHist[i]->Write();	
	}
		
	// Saving everything...
	output->cd();
	hist4MuMUnconstrained->Write();
	hist4ElMUnconstrained->Write();
	hist2El2MuMUnconstrained->Write();
	hist2Mu2ElMUnconstrained->Write();
	
	hist4MuMFSR->Write();
	hist4ElMFSR->Write();
	hist2El2MuMFSR->Write();
	hist2Mu2ElMFSR->Write();
	
	hist4MuMConstrained->Write();
	hist4ElMConstrained->Write();
	hist2El2MuMConstrained->Write();
	hist2Mu2ElMConstrained->Write();
	muInteractionHist->Write();
	muInteractionHistPileupW->Write();

	countingHist->Write();
	// creating a folder, which is the name of the dataContainer
	output->cd(); 
	output->mkdir(fileNamePart);

	output->Write();

	output->Close();

}

