#ifndef CHARGEDLEPTON_H
#define CHARGEDLEPTON_H

#include <D3PDReader/MuonD3PDObject.h>
#include <D3PDReader/ElectronD3PDObject.h>
#include <D3PDReader/JetD3PDObject.h>

#include <MyAnalysis/macroDef.h>

#include "ZMassConstraint/CovMatrixTools.h"
#include "egammaAnalysisUtils/EnergyRescaler.h"
#include "egammaAnalysisUtils/EnergyRescalerUpgrade.h"
#include "ElectronPhotonFourMomentumCorrection/egammaEnergyCorrectionTool.h"

#include <TTree.h>
#include <TString.h>
#include <TLorentzVector.h>

using namespace CLHEP;

class ChargedLepton 
{
	public :
		// Main Variables
		Int_t flavor; //Muon or Electron
		Int_t type; // What type of lepton is it
		Int_t charge; // Charge of the lepton for easier excess afterwards.
		Int_t index;

		Double_t ptCone20; // value in the ptCone20 var
		Double_t ptCone20Correction; // value to subtract from the ptCone20 if conditions are meet. See HiggsAnalysis::getTrackPT
	
		Double_t smearVal;
		Double_t lepEff;
		Double_t lepMass;

		// Lepton ID
		Int_t lepID;

		// Momentum Error Term
		Double_t covMomErr;

		// Resolution Term for electron
		Double_t resolutionEl;

		// Cluster PT for electron
		Double_t clPt;

		// For william's systematics
		Float_t id_pt_unsmeared;
		Float_t me_pt_unsmeared;
		Float_t cb_pt_unsmeared;
		// For Graham's systematics
		Float_t E_sys[doSys::Nom + 1];	
		Float_t E_nom;

		// For categorization
		Int_t lepType;
		Int_t truthParentType;

		//Main TLorentz Vector
		TLorentzVector *m_momentum_main;
		TLorentzVector m_momentum; // For FSR and ZMass Constraint... I messed up with this
		TLorentzVector m_momentumMS; // For MS muons
		TLorentzVector m_momentumID; // For ID muons
		TLorentzVector m_momentumBDT; // BDT tool need a TVL that is GeV
		TLorentzVector m_momentumTrig; // For Trigger Efficiency Calculation
		TLorentzVector m_momentumTruthRecoBare; // For reco Truth bare calculations	
		TLorentzVector m_momentumTruthRecoBorn; // For reco Truth born calculations	
		TLorentzVector m_momentumTruthTrueBorn; // For true turht born calculations

		// Error Matrix
		TMatrixD errorMatrix;
		HepMatrix errorMatrixHep;	
		TMatrixD errorMatrixMS;
		HepMatrix errorMatrixHepMS;
		TMatrixD errorMatrixID;
		HepMatrix errorMatrixHepID;
		// Rescaler
		AtlasRoot::egammaEnergyCorrectionTool *elRescale;

		D3PDReader::MuonD3PDObjectElement *mu;
		D3PDReader::ElectronD3PDObjectElement *el;
		D3PDReader::JetD3PDObjectElement *jets;

		// Getter
		TLorentzVector *get4Momentum() {return m_momentum_main;}
		TLorentzVector get4MomentumNoP(Int_t muType = muonType::CB) {
			if(muType == muonType::CB) return m_momentum;
			if(muType == muonType::MS) return m_momentumMS;
			if(muType == muonType::ID) return m_momentumID;
			else return m_momentum;
		}
		TLorentzVector get4MomentumBDT() {return m_momentumBDT;}
		TLorentzVector get4MomentumTrigSF() {return m_momentumTrig;}				
		
		Int_t getFlavor() {return flavor;}		
		Int_t getCharge() {return charge;}		
		Int_t getType() {return type;}		
		
		Double_t getPTCone20() {return ptCone20;}	
		Double_t getPRCone20Corr() {return ptCone20Correction;}		
		
		TMatrixD getCovMatrix(Int_t muType = muonType::CB)
		{
			if(muType == muonType::CB) return errorMatrix;
			if(muType == muonType::MS) return errorMatrixMS;
			if(muType == muonType::ID) return errorMatrixID;
			else return errorMatrix;
		}		
		HepMatrix getHepCovMatrix(Int_t muType = muonType::CB) {
			if(muType == muonType::CB) return errorMatrixHep;
			if(muType == muonType::MS) return errorMatrixHepMS;
			if(muType == muonType::ID) return errorMatrixHepID;
			else return errorMatrixHep;;
		}		
		
		Double_t getLepEff(){return lepEff;}

<<<<<<< HEAD
		// For Filling the Covariance Matrix	
=======
		// For Filling the Covariance Martix	
>>>>>>> bc7b9ddaf72f0a41dfe1bb5d9068cc4b03444c0d
		void SetElRescale(AtlasRoot::egammaEnergyCorrectionTool *telRescale){elRescale = telRescale;}

		void FillCovMatrix(Int_t runNumber_sf);
		// Constructor & Destructor
		
		// Muon Constructor Differentiates between the staco CB and standalone 
		ChargedLepton(D3PDReader::MuonD3PDObjectElement *tMu, Int_t tType, Double_t tLepEff, Int_t tindex, Bool_t isMC);
		ChargedLepton(D3PDReader::ElectronD3PDObjectElement *tEl, Int_t tType, Double_t tLepEff, Double_t tSmearVal, Int_t tindex, Bool_t isMC, Double_t tresolution = -1, Double_t clPt = -1);
		ChargedLepton(D3PDReader::JetD3PDObjectElement *tJets, Int_t tType, Int_t tindex, Bool_t isMC);
		~ChargedLepton();

		D3PDReader::MuonD3PDObjectElement* GetMuon() {return mu;}
		D3PDReader::ElectronD3PDObjectElement* GetElectron (){return el;}
		D3PDReader::JetD3PDObjectElement* GetJets (){return jets;}

		
	private :
};

#endif



