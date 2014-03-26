#ifndef CORRECTION_H
#define CORRECTION_H

#include "MyAnalysis/macroDef.h"

#include "D3PDReader/Event.h"
#include "D3PDReader/MuonD3PDObject.h"
#include "D3PDReader/ElectronD3PDObject.h"
#include "D3PDReader/PhotonD3PDObject.h"
#include "D3PDReader/JetD3PDObject.h"

#include "MuonMomentumCorrections/SmearingClass.h"
#include "MuonMomentumCorrections/MuonResolutionAndMomentumScaleFactors.h"
#include "egammaAnalysisUtils/EnergyRescaler.h"
#include "egammaAnalysisUtils/EnergyRescalerUpgrade.h"
#include "ApplyJetCalibration/ApplyJetCalibration.h"

#include "ElectronEfficiencyCorrection/TElectronEfficiencyCorrectionTool.h"
#include "MuonEfficiencyCorrections/AnalysisMuonConfigurableScaleFactors.h"
#include "MuonEfficiencyCorrections/AnalysisMuonEfficiencyScaleFactors.h"

#include "PileupReweighting/TPileupReweighting.h"

#include "ElectronPhotonFourMomentumCorrection/egammaEnergyCorrectionTool.h"
#include "egammaFourMomentumError/egammaFourMomentumError.h"

#include "MyAnalysis/ChargedLepton.h"
#include "MyAnalysis/macroDef.h"

#include <TH2F.h>
#include <TFile.h>
#include <TAxis.h>
#include <TRandom3.h>
#include <TLorentzVector.h>

#include "math.h"
#include "vector"

class Correction 
{
	public:
		// Constructor & Destructor
		Correction(Int_t tdataYear, Bool_t tIsMC, Bool_t tSmearMC, Bool_t tScaleData, 
				Bool_t tScaleCrack, Bool_t tJetCalibration, Bool_t tScaleMuon, Bool_t tscaleEfficiency, Bool_t tuseMoriond, 
				Int_t tcurrCollection, Int_t tcurDataCollection, Bool_t tepCombination, Root::TPileupReweighting* tpileupTool);
		~Correction();
		// Vars to Control Smearing
		Int_t dataYear;
		Int_t currCollection;
		Int_t currDataCollection;
		Bool_t isMC;
		Bool_t smearMC;
		Bool_t scaleData;
		Bool_t scaleCrack;
		Bool_t jetCalibration;
		Bool_t scaleMuon;
		Bool_t scaleEfficiency;
		Bool_t epCombination;
		Bool_t useMoriond;

		Bool_t isDebugCall;

		// To Print
		Int_t printInfoMuSmear;
		Int_t printInfoElSmear;		
		Int_t printInfoPhSmear;		
		Int_t printInfoD0Z0Smear;
		Int_t printInfoJetCal;
		Int_t printInfoElEff;
		Int_t printInfoMuEff;

		// Before EP combination cl_pt for FSR search
		vector <float> bfEP_cl_pt;
		vector <float> bfEP_cl_Et;

		// D0/Z0 smear
		void InitD0Z0Smear();
		void SmearD0Z0(D3PDReader::ElectronD3PDObject * el, Int_t EventNumber);
		void SmearD0Z0(D3PDReader::MuonD3PDObject * mu, Int_t EventNumber, Int_t type);

		void debugCall(){isDebugCall = true;}

		// lepton and jets Correction
		Root::TPileupReweighting* pileupTool;
		MuonSmear::SmearingClass* muSmear;
		Analysis::MuonResolutionAndMomentumScaleFactors *mResolSF;
  		//egRescaler::EnergyRescalerUpgrade *elRescale;
  		AtlasRoot::egammaEnergyCorrectionTool* elRescale;		
  		JetAnalysisCalib::JetCalibrationTool *calibrationJES;	
 		// Efficiency Calculations
		Root::TElectronEfficiencyCorrectionTool *egSFClassID;
  		Root::TElectronEfficiencyCorrectionTool *egSFClassReco;
		// Staco Efficiency Calculations
		Analysis::AnalysisMuonConfigurableScaleFactors* StacoSCF;
		// Staco SA Efficiency Calculations
		Analysis::AnalysisMuonConfigurableScaleFactors* StacoSASCF;
		// CaloMuon Efficiency Calculations
		Analysis::AnalysisMuonConfigurableScaleFactors* CaloMuSCF;

		// EP combination
		egammaFourMomentumError *EPcombination;
		
		void InitMuonSmear();
		void InitElectronSmear(Int_t electronCollection);	
		void InitJetCal();				
		void SmearMuon(D3PDReader::MuonD3PDObject * mu, Int_t EventNumber, Int_t type);
		void CalibrateJet(D3PDReader::JetD3PDObject * jet, Int_t dataYear, Double_t rhoKt4EM, Double_t mu, Double_t NPV);
		vector<Double_t> SmearElectron(D3PDReader::ElectronD3PDObject * el, Int_t EventNumber, Int_t RunNumber, Int_t sysVar = doSys::Nom);
		vector<Double_t> SmearPhoton(D3PDReader::PhotonD3PDObject * ph, Int_t EventNumber, Int_t RunNumber);

		void ClearVars();

		vector<Double_t> electronEff;
		vector<Double_t> electronEpErr;
		vector<Double_t> muonCaloEff;
		vector<Double_t> muonStacoEff;

		AtlasRoot::egammaEnergyCorrectionTool* getElRescale(){return elRescale;}
		
	private:
		  // d0 smearing
  		TH2F *smearD0[3];
		TRandom3 smearD0Rand;
		TAxis *smearD0X;
		TAxis *smearD0Y;
		
		// z0 smearing
		TH2F *smearZ0[3];
		TRandom3 smearZ0Rand;
		TAxis *smearZ0X;
		TAxis *smearZ0Y;
		Double_t GetD0SmearSigma(Int_t EventNumber, Int_t nBL, Double_t pt, Double_t eta, Int_t index);
		Double_t GetZ0SmearSigma(Int_t EventNumber, Int_t nBL, Double_t pt, Double_t eta, Int_t index);
		
};

#endif



