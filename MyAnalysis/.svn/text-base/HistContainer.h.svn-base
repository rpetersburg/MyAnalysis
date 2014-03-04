#ifndef HISTCONTAINER_H
#define HISTCONTAINER_H


#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TString.h>

class HistContainer 
{
	public :
		// weight for normalization
		Double_t weight;
		// Histograms
		// Plotting Varibles
		TH1D *hist4MuMUnconstrained;
		TH1D *hist4ElMUnconstrained;
		TH1D *hist2El2MuMUnconstrained;
		TH1D *hist2Mu2ElMUnconstrained;
		TH1D *hist4MuMFSR;
		TH1D *hist4ElMFSR;
		TH1D *hist2El2MuMFSR;
		TH1D *hist2Mu2ElMFSR;
		TH1D *hist4MuMConstrained;
		TH1D *hist4ElMConstrained;
		TH1D *hist2El2MuMConstrained;
		TH1D *hist2Mu2ElMConstrained;	

		// For Counting Histrogram
		TH1D *cutPassHist;
		TH1D *cutMuPassHist;
		TH1D *cutElPassHist;
		TH1D *cutJetsPassHist;
		TH1D *cut4MuPassHist;
		TH1D *cut4ElPassHist;
		TH1D *cut2L2LPassHist;
		TH1D *truthQuadHist;
		// for Weights
		TH1D *cutPassHistW;
		TH1D *cut4MuPassHistW;
		TH1D *cut4ElPassHistW;
		TH1D *cut2L2LPassHistW;
		TH1D *cutlleePassHistW;

		// Information Storing Function For the cuts
		// Muons
		TH1D **muAuthorHist;
		TH1D **muPTHist;
		TH1D **muEtaHist;
		TH1D **muD0Hist;
		TH1D **muZ0Hist;
		TH1D **muOverlapHist;
		// Electrons
		TH1D *elAuthorHist;
		TH1D *elETHist;
		TH1D *elEtaHist;
		TH1D *elZ0Hist;
		// QuadLepton
		TH1D ***quadInitPtHist;
		TH1D ***quadFinalPtHist;
		TH1D **Z1MassHist;
		TH2D **Z2M4lHist;
		TH1D **JPsiHist;
		TH1D **TrackIsoHist;
		TH1D **muCaloIsoHist;
		TH1D **elCaloIsoHist;
		TH1D **muSigHist;
		TH1D **elSigHist;

		// General
		TH1D *muInteractionHist; 		
		TH1D *muInteractionHistPileupW; 		

		// Constructor & Destructor
		HistContainer(Int_t nCut, Int_t nMuCut, Int_t nElCut, Int_t nJetsCut, Int_t nCH);
		~HistContainer();

		void SaveHist(TFile *output, TString fileNamePart, TH1F* countingHist);
		
	private :
		
};

#endif



