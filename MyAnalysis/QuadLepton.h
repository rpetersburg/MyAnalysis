#ifndef QUADLEPTON_H
#define QUADLEPTON_H


#include <MyAnalysis/macroDef.h>
#include <MyAnalysis/ChargedLepton.h>
#include <MyAnalysis/DiLepton.h>

#include <TLorentzVector.h>

#include "egammaAnalysisUtils/EnergyRescaler.h"
#include "egammaAnalysisUtils/EnergyRescalerUpgrade.h"

#include "ElectronPhotonFourMomentumCorrection/egammaEnergyCorrectionTool.h"

#include <vector>
using namespace CLHEP;

class QuadLepton 
{
	public :
		// Main Variables
		Int_t quadType; 
		Bool_t neutral; 
		Int_t prodChannel;
		Int_t quadTypeBR;

		//Main TLorentz Vector
		TLorentzVector *m_momentum_main;
		TLorentzVector m_momentum;		
		TLorentzVector fsr_momentum;		

		// Vector to have easier access to all the leptons
		vector<ChargedLepton *> leptonInfo;
		
		vector<TLorentzVector> leptonLorentz;
		vector<TMatrixD> leptonCovMatrix;
		vector<HepMatrix> leptonHepCovMatrix;	
		// For FSR
		vector<TLorentzVector> leptonFSRLorentz;
		vector<TMatrixD> leptonFSRCovMatrix;
		vector<HepMatrix> leptonHepFSRCovMatrix;
		// For ZMass
		vector<TLorentzVector> leptonZMassFSRLorentz;
		vector<TMatrixD> leptonZMassCovMatrix;

		vector<TLorentzVector> leptonLorentzID;
		vector<TMatrixD> leptonCovMatrixID;
		vector<HepMatrix> leptonHepCovMatrixID;
		vector<TLorentzVector> leptonLorentzMS;
		vector<TMatrixD> leptonCovMatrixMS;
		vector<HepMatrix> leptonHepCovMatrixMS;

		// For FSR
		vector<TLorentzVector> leptonFSRLorentzID;
		vector<TMatrixD> leptonFSRCovMatrixID;
		vector<HepMatrix> leptonHepFSRCovMatrixID;
		vector<TLorentzVector> leptonFSRLorentzMS;
		vector<TMatrixD> leptonFSRCovMatrixMS;
		vector<HepMatrix> leptonHepFSRCovMatrixMS;

		// For ZMass
		vector<TLorentzVector> leptonZMassFSRLorentzID;
		vector<TMatrixD> leptonZMassCovMatrixID;
		vector<TLorentzVector> leptonZMassFSRLorentzMS;
		vector<TMatrixD> leptonZMassCovMatrixMS;



		// For CR
		vector<Bool_t> looseElectron;
		vector<Bool_t> trackIso;
		vector<Bool_t> caloIso;
		vector<Bool_t> d0Sig;
		vector<Double_t> valTrackIso;
		vector<Double_t> valCaloIso;
		vector<Double_t> valD0Sig;

		// Mass information
		Double_t mass;
		Double_t massErr;
		// ID
		Double_t massID;
		Double_t massErrID;
		// MS
		Double_t massMS;
		Double_t massErrMS;

		Double_t massFSR;
		Double_t massErrFSR;
		Double_t massZ1FSR;
		Double_t massZ2FSR;
		Int_t 	 fsrType;

		Double_t massZMassCons;
		Double_t massErrZmassCons;
		Double_t massZ1ZMassCons;
		Double_t massZ2ZMassCons;

		//ID
		Double_t massZMassConsID;
		Double_t massErrZmassConsID;
		//MS
		Double_t massZMassConsMS;
		Double_t massErrZmassConsMS;


		TLorentzVector sum_unconstrained;
		TLorentzVector sum_fsr;
		TLorentzVector sum_constrained;
		TLorentzVector truthVec;

		//Production Angles
		Float_t cthstr;
		Float_t phi1;
		Float_t cth1;
		Float_t cth2;
		Float_t phi;

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

		// Truth information
		Double_t m4l_truth_born; 
		Double_t mZ1_truth_born;
		Double_t mZ2_truth_born;
		Double_t m4l_truth_bare; 
		Double_t mZ1_truth_bare; 
		Double_t mZ2_truth_bare;
		Double_t m4l_truth_dressed; 
		Double_t mZ1_truth_dressed; 
		Double_t mZ2_truth_dressed;

		// For Jets
		ChargedLepton *leadingJet;
		ChargedLepton *subLeadingJet;
		ChargedLepton *thirdJet;
		Int_t	n_jets;
		Float_t dijet_invmass;
		Float_t dijet_deltaeta;
		Float_t leading_jet_pt;
		Float_t leading_jet_eta;
		Float_t subleading_jet_pt;
		// Truth
		Int_t 	n_jets_truth_bare;
		Float_t	leading_jet_pt_truth_bare;

		// Truth
		Int_t 	n_jets_fid;
		Float_t	leading_jet_pt_fid;
		Int_t 	n_jets_truth_fid;
		Float_t	leading_jet_pt_truth_fid;

		Float_t BDT_discriminant_VBF;
		Float_t BDT_discriminant_HadVH;

		// BDT
		Float_t KD_discriminant;
		Float_t BDT_discriminant;
		Float_t BDTGuass_discriminant;
		Float_t ptSysupFac;
		Float_t ptSysdownFac;

		// NPV
		Int_t npv;	

		// BCHCut
		Int_t BCHCutMedium;
		Int_t BCHCutTight;
		
		// Type of calibariton
		Int_t calib;
		
		// Flag for overlap in data stream
		Int_t flagQuad;

		// For categories
		ChargedLepton* leadingExtraLep;
		ChargedLepton* subleadingExtraLep;
		// Functions
		DiLepton *Z1;
		DiLepton *Z2;

		// For photon covariance Matrix
		AtlasRoot::egammaEnergyCorrectionTool *elRescale;
		// Getter
		TLorentzVector *get4Momentum() {return m_momentum_main;}
		TLorentzVector get4MomentumNoP() {return m_momentum;}

		DiLepton *getZ1() {return Z1;}
		DiLepton *getZ2() {return Z2;}
		Bool_t IsNeutral() {return neutral;}
		Int_t getQuadType() {return quadType;}
		Int_t getQuadTypeBR() {return quadTypeBR;}
		
		vector<ChargedLepton *> getLepton() {return leptonInfo;}
		
		vector<TLorentzVector> getLeptonLorentz(Int_t muType = muonType::CB) {
			if(muType == muonType::CB) return leptonLorentz;
			if(muType == muonType::MS) return leptonLorentzMS;
			if(muType == muonType::ID) return leptonLorentzID;
			else return leptonLorentz;
		}
		vector<TMatrixD> getLeptonCovMatrix(Int_t muType = muonType::CB) {
			if(muType == muonType::CB) return leptonCovMatrix;
			if(muType == muonType::MS) return leptonCovMatrixMS;
			if(muType == muonType::ID) return leptonCovMatrixID;
			else return leptonCovMatrix;
		}
		// For FSR
		vector<TLorentzVector> getLeptonFSRLorentz() {return leptonFSRLorentz;}
		vector<TMatrixD> getLeptonFSRCovMatrix() {return leptonFSRCovMatrix;}
		// For Z mass
		vector<TLorentzVector> getLeptonZMassLorentz() {return leptonZMassFSRLorentz;}
		vector<TMatrixD> getLeptonZMassCovMatrix() {return leptonZMassCovMatrix;}


		Double_t getMass() {return mass;}
		Double_t getMassErr() {return massErr;}
		
		Double_t getMassFSR() {return massFSR;}
		Double_t getMassErrFSR() {return massErrFSR;}
		Double_t getZ1MassFSR() {return massZ1FSR;}
		Double_t getZ2MassFSR() {return massZ2FSR;}
		
		Double_t getMassZMassCons() {return massZMassCons;}
		Double_t getMassErrZMassCons() {return massErrZmassCons;}
		Double_t getZ1MassZMassCons() {return massZ1ZMassCons;}
		Double_t getZ2MassZMassCons() {return massZ2ZMassCons;}
	
		Int_t getProductionChannel() {return prodChannel;}

		TLorentzVector getFSR4MomentumNoP() {return fsr_momentum;}

		// setters
		void SetMass(Double_t tmass) {mass = tmass;}
		void SetMassErr(Double_t tmassErr) {massErr = tmassErr;}
		
		void SetMassFSR(Double_t tmassFSR) {massFSR = tmassFSR;}
		void SetMassErrFSR(Double_t tmassErrFSR) {massErrFSR = tmassErrFSR;}
		void SetZ1FSRMass(Double_t tmass) {massZ1FSR = tmass;}	
		void SetZ2FSRMass(Double_t tmass) {massZ2FSR = tmass;}	
		
		void SetMassZMass(Double_t tmassZMassCons, Int_t muType) 
		{
			if(muType == muonType::CB) massZMassCons = tmassZMassCons;
			if(muType == muonType::MS) massZMassConsMS = tmassZMassCons;
			if(muType == muonType::ID) massZMassConsID = tmassZMassCons;
		}
		void SetMassErrZMass(Double_t tmassErrZmassCons, Int_t muType) 
		{
			if(muType == muonType::CB) massErrZmassCons = tmassErrZmassCons;
			if(muType == muonType::MS) massErrZmassConsMS = tmassErrZmassCons;
			if(muType == muonType::ID) massErrZmassConsID = tmassErrZmassCons;
		}
		void SetZ1MassZMass(Double_t tmassZ1ZMassCons, Int_t muType) 
		{
			if(muType == muonType::CB) massZ1ZMassCons = tmassZ1ZMassCons;
		}
		void SetZ2MassZMass(Double_t tmassZ2ZMassCons, Int_t muType) 
		{
			if(muType == muonType::CB)massZ2ZMassCons = tmassZ2ZMassCons;
		}
		
		void SetZmassLorentz(vector<TLorentzVector> tleptonZMassFSRLorentz) {leptonZMassFSRLorentz = tleptonZMassFSRLorentz;}
		void SetZmassCovMat(vector<TMatrixD> tleptonZMassCovMatrix) {leptonZMassCovMatrix = tleptonZMassCovMatrix;}
		
		void SetProductionChannel (Int_t tprodChannel) {prodChannel = tprodChannel;}
		void FillFSRCorr(TLorentzVector momFSR, Int_t runNumber_sf, Bool_t isZ1, Bool_t isZ2, PATCore::ParticleType::Type type = PATCore::ParticleType::Electron);

		// For Filling the Covariance Martix			
		void SetElRescale(AtlasRoot::egammaEnergyCorrectionTool *telRescale);
		void FillCovMatrix(Int_t runNumber_sf);


		// Constructor & Destructor
		QuadLepton();
		QuadLepton(DiLepton *ZOne, DiLepton *ZTwo);
		~QuadLepton();

		void Set(DiLepton *ZOne, DiLepton *ZTwo);

		
	private :
		void Reset();
};

#endif



