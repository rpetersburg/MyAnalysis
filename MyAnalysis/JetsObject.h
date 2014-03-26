#ifndef JETSOBJECT_H
#define JETSOBJECT_H

#include <D3PDReader/JetD3PDObject.h>
#include <MyAnalysis/macroDef.h>
#include <MyAnalysis/ChargedLepton.h>

#include <TTree.h>
#include <TString.h>
#include <TMath.h>
#include <TLorentzVector.h>

#include <math.h>
#include <vector>

class JetsObject 
{
	public :

		Int_t dataYear;
		Int_t eventNumber;

		// Storage of Variables
		vector<ChargedLepton *> jetsInitEvent; // Initial Storage
		vector<ChargedLepton *> jetsAuthor;	// Before Overlap but after all other cuts			
		vector<ChargedLepton *> jetsBfOverlap;	// Before Overlap but after all other cuts	
		vector<ChargedLepton *> jetsOverlapGoodEvent; // Events that have Passed the Jets cuts plus overlap
		vector<ChargedLepton *> jetsEvent;

		// Constructor & Destructor
		
		// Jets Constructor Differentiates between the staco CB and standalone 
		JetsObject(Int_t year);

		// Filling in the vector
		void FillJets(D3PDReader::JetD3PDObject * jets_branch, Int_t type, Bool_t isMC, Int_t teventNumber);
		void FillJetsTruth(D3PDReader::JetD3PDObject * jets_branch, Int_t type, Bool_t isMC);
		
		Bool_t JetsCut(Int_t *cutJetsPass, Int_t run);
		
		vector<ChargedLepton *> getJetsVec() {return jetsBfOverlap;}
		void clearVars();

	private:
		Bool_t CutAntiKt4TopoEM(D3PDReader::JetD3PDObjectElement *jets_i, Int_t *cutJetsPass, Int_t run);
		Bool_t CutAntiKt4TopoEMTruth(D3PDReader::JetD3PDObjectElement *jets_i);
		Bool_t CutAntiKt4TopoEM_Fid(D3PDReader::JetD3PDObjectElement *jets_i, Int_t run);
		Bool_t CutAntiKt4TopoEMTruth_Fid(D3PDReader::JetD3PDObjectElement *jets_i);
		
		Double_t DeltaR (Double_t eta_1, Double_t phi_1, Double_t eta_2, Double_t phi_2);
			


};

#endif



