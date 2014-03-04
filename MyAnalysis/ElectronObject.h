#ifndef ELECTRONOBJECT_H
#define ELECTRONOBJECT_H

#include "D3PDReader/ElectronD3PDObject.h"
#include "MyAnalysis/macroDef.h"
#include "MyAnalysis/ChargedLepton.h"
#include "MyAnalysis/HistContainer.h"

//#include "egammaAnalysisUtils/H4l2011Defs.h"
//#include "egammaAnalysisUtils/MultiLeptonMenu.h"

#include "ElectronPhotonSelectorTools/TElectronLikelihoodTool.h"
#include "ElectronPhotonSelectorTools/TElectronIsEMSelector.h"
#include "ElectronPhotonSelectorTools/TElectronMultiLeptonSelector.h"

#include <TTree.h>
#include <TH1F.h>
#include <TString.h>
#include <TMath.h>
#include <TString.h>
#include "TPython.h"

#include <math.h>
#include <vector>

class ElectronObject 
{
	public :

		// For Cut Flow
		Int_t dataYear;
		Bool_t useLikelihood;
  		//H4l2011Defs* elID2011;
  		Root::TElectronIsEMSelector* elID2011new;
 		Root::TElectronLikelihoodTool* elID2012;		
	  	Root::TElectronMultiLeptonSelector* ml_2013;
	
  		//MultiLeptonMenu* ml_2013;
		
		// Storage of Variables
		vector<ChargedLepton *> elInitEvent; // Initial Storage
		vector<ChargedLepton *> elBfOverlap; // Events that have Passed the Electron but not overlap		
		vector<ChargedLepton *> elBfCCOverlap; // Events that have Passed the Electron cuts plus ee overlap	but not eecc
		vector<ChargedLepton *> elOverlapGoodEvent; // Events that have Passed the EKectron cuts plus overlap
		vector<ChargedLepton *> elEvent;

		// To keep the calibration
		vector<Double_t> elSmearVal;		

		void clearVars();

		// Constructor & Destructor
		ElectronObject(Int_t year, Bool_t tuseLikelihood);

		// Filling in the vector
		void FillElectron(D3PDReader::ElectronD3PDObject * el_branch, Int_t type, vector<Double_t> telSmearVal, vector<Double_t> elEff, vector<Double_t> elRes, vector<Double_t> elClPT, Bool_t isMC);
		Bool_t ElectronCut(Int_t *cutElPass, Int_t Npv, Bool_t useRelaxedLoose = false);

		vector<ChargedLepton *> getElectronVec() {return elOverlapGoodEvent;}

		// Histograms for Storage
		HistContainer *Hist;		
	
		void SetHist(HistContainer *curr_Hist);
		
		Bool_t AuthorCut(D3PDReader::ElectronD3PDObjectElement *el_i, Double_t smearVal_i, Int_t ip, Bool_t useRelaxed = false);

		
	private:
		// Overlap
		void RemoveOverlap(Int_t * cutElPass, Int_t NPV);
		Bool_t OverlapEE(D3PDReader::ElectronD3PDObjectElement *el_curr, Int_t NPV);
		Bool_t OverlapEECC(D3PDReader::ElectronD3PDObjectElement *el_curr, Int_t NPV);
		

		// Helpers
		Double_t DeltaR (Double_t eta_1, Double_t phi_1, Double_t eta_2, Double_t phi_2);
			
		// Specific Cuts
		Bool_t CutGSF(D3PDReader::ElectronD3PDObjectElement *el_i, Int_t *cutElPass, Double_t z0Cut, Double_t smearVal_i, Int_t Npv, Bool_t useRelaxed = false);
	
};

#endif



