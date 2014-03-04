#ifndef DILEPTON_H
#define DILEPTON_H


#include <MyAnalysis/macroDef.h>
#include <MyAnalysis/ChargedLepton.h>
#include <TLorentzVector.h>
#include "ElectronPhotonFourMomentumCorrection/egammaEnergyCorrectionTool.h"

using namespace CLHEP;

class DiLepton 
{
	public :
		// Main Variables
		Int_t flavor; //Muon or Electron
		Bool_t neutral; // Gets the overall charged of the two leptons
		Bool_t hasFSR;
		Int_t type;

		//Main TLorentz Vector
		TLorentzVector *m_momentum_main;
		TLorentzVector m_momentum;
		TLorentzVector FSR_momentum;
		TMatrixD FSRError;
		HepMatrix HepFSRError;

		ChargedLepton *lepPlus;
		ChargedLepton *lepNeg;

		// Vars for Dilep analysis
		// Mass information
		Double_t mass;
		Double_t massErr;
		Double_t massZMassCons;
		Double_t massErrZmassCons;
		
		// Vector to have easier access to all the leptons
		vector<ChargedLepton *> leptonInfo;
		vector<TLorentzVector> leptonLorentz;		
		
		vector<TMatrixD> leptonCovMatrix;
		vector<HepMatrix> leptonHepCovMatrix;
		vector<TLorentzVector> getLeptonLorentz() {return leptonLorentz;}
		vector<TMatrixD> getLeptonCovMatrix() {return leptonCovMatrix;}
		// For ZMass
		vector<TLorentzVector> leptonZMassFSRLorentz;
		vector<TMatrixD> leptonZMassCovMatrix;
		void SetZmassLorentz(vector<TLorentzVector> tleptonZMassFSRLorentz) {leptonZMassFSRLorentz = tleptonZMassFSRLorentz;}
		void SetZmassCovMat(vector<TMatrixD> tleptonZMassCovMatrix) {leptonZMassCovMatrix = tleptonZMassCovMatrix;}
		vector<TLorentzVector> getLeptonZMassLorentz() {return leptonZMassFSRLorentz;}
		vector<TMatrixD> getLeptonZMassCovMatrix() {return leptonZMassCovMatrix;}

		// Getter
		TLorentzVector *get4Momentum() {return m_momentum_main;}
		TLorentzVector get4MomentumNoP() {return m_momentum;}
		
		// For photon covariance Matrix
		AtlasRoot::egammaEnergyCorrectionTool *elRescale;
		
		ChargedLepton *getLepPlus() {return lepPlus;}
		ChargedLepton *getLepNeg() {return lepNeg;}
		Bool_t IsNeutral() {return neutral;}
		Int_t getFlavor() {return flavor;}

		Bool_t getHasFSR() {return hasFSR;}
		TLorentzVector getFSR4Momentum () {return FSR_momentum;}
		TMatrixD getFSRError() {return FSRError;}
		HepMatrix getHepFSRError() {return HepFSRError;}

		vector<ChargedLepton *> getLepton() {return leptonInfo;}



		// vars
		// Weights
		Double_t weight;
		Double_t weight_corr;
		Double_t weight_lumi;
		Double_t zVertWeight;
		Double_t ggFWeight;		
		Double_t JHUWeight;		
		Double_t pileupWeight;
		Double_t ttBarVeto;
		Double_t zzOverlapVeto;
		Double_t zBBVeto;
		Double_t trigEff;
		Double_t lepEff;
		Double_t mcEventWeight;
		Double_t crossSection;
		Double_t branchRatio;
		Double_t lumi;

		// Setters
		void setHasFSR(Bool_t tHasFSR){hasFSR = tHasFSR;}
		void setFSR4Momentum (TLorentzVector tFSR_momentum){FSR_momentum = tFSR_momentum;}
		void setFSRError (TMatrixD tFSRError){FSRError.ResizeTo(tFSRError); FSRError = tFSRError;}
		void setHepFSRError (HepMatrix tHepFSRError){ HepFSRError = tHepFSRError;}
		
		// For Filling the Covariance Martix			
		void SetElRescale(AtlasRoot::egammaEnergyCorrectionTool *telRescale);
		void FillCovMatrix(Int_t runNumber_sf);

		// Constructor & Destructor
		DiLepton();
		~DiLepton();		
		DiLepton(ChargedLepton *lepOne, ChargedLepton *lepTwo);

		void Set(ChargedLepton *lepOne, ChargedLepton *lepTwo);

		Bool_t IsOverlap (DiLepton *toCompare);
		
	private :
		void Reset();
};

#endif



