#ifndef MUONOBJECT_H
#define MUONOBJECT_H

#include <D3PDReader/MuonD3PDObject.h>
#include <MyAnalysis/macroDef.h>
#include <MyAnalysis/ChargedLepton.h>
#include <MyAnalysis/HistContainer.h>

#include <TTree.h>
#include <TH1F.h>
#include <TString.h>
#include <TMath.h>

#include <math.h>
#include <vector>

class MuonObject 
{
	public :

		Int_t dataYear;
		Int_t cbAuthor;
		Int_t saAuthor;
		Int_t caloAuthor;

		Bool_t looseCutZ4l;

		// Storage of Variables
		vector<ChargedLepton *> muInitEvent; // Initial Storage
		vector<ChargedLepton *> muBfOverlap;	// Before Overlap but after all other cuts	
		vector<ChargedLepton *> muOverlapGoodEvent; // Events that have Passed the Muon cuts plus overlap
		vector<ChargedLepton *> muEvent;

		// Constructor & Destructor
		// Muon Constructor Differentiates between the staco CB and standalone 
		MuonObject(Int_t year, Bool_t tLooseCutZ4l = false);
		~MuonObject();
		

		// Filling in the vector
		void FillMuon(D3PDReader::MuonD3PDObject * mu_branch, Int_t type, vector<Double_t> muEff, Bool_t isMC);
		Bool_t MuonCut(Int_t *cutMuPass);
		
		vector<ChargedLepton *> getMuonVec() {return muOverlapGoodEvent;}

		void SetHist (HistContainer *curr_Hist);
		void clearVars();		

		// Histrograms
		HistContainer *Hist;
		
	private:
		// Helpers
		void RemoveOverlap(Int_t *cutMuPass);
		Bool_t OverlapCalo(D3PDReader::MuonD3PDObjectElement *mu_curr);
		Bool_t OverlapStandAlone(D3PDReader::MuonD3PDObjectElement *mu_curr);


		Double_t DeltaR (Double_t eta_1, Double_t phi_1, Double_t eta_2, Double_t phi_2);
			
		// Specific Cuts
		Bool_t CutStaco(D3PDReader::MuonD3PDObjectElement *mu_i, Int_t *cutMuPass, Double_t d0Cut, Double_t z0Cut);
		Bool_t CutCalo(D3PDReader::MuonD3PDObjectElement *mu_i, Int_t *cutMuPass, Double_t d0Cut, Double_t z0Cut);
		Bool_t CutStandAlone(D3PDReader::MuonD3PDObjectElement *mu_i, Int_t *cutMuPass, Double_t d0Cut, Double_t z0Cut);

};

#endif



