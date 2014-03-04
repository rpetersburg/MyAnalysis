#ifndef OUTPUTTREESYS_H
#define OUTPUTTREESYS_H

#include <TTree.h>
#include <TString.h>
#include <TLorentzVector.h>

#include "D3PDReader/Event.h"

#include "MyAnalysis/QuadLepton.h"
#include "MyAnalysis/macroDef.h"

#include <stdlib.h>
#include <math.h>
#include <iomanip>
#include <string>
#include <TH1F.h>


class OutputTreeSys 
{
	public :
		OutputTreeSys();

		~OutputTreeSys();
		
		void clearVars();
		
		void fillTree(D3PDReader::Event *event, QuadLepton * higgs, Int_t type, Bool_t isMC);
		void saveTrees(TString filePath, TH1F* countingHist, TString sampleName);
			
		void bookTree(TTree *tree);
	private :
		// TTree for Each channel
		TTree * mu4Tree;
		TTree * el4Tree;
		TTree * mu2el2Tree;
		TTree * el2mu2Tree;

		TTree * mu4ggFTree;
		TTree * el4ggFTree;
		TTree * mu2el2ggFTree;
		TTree * el2mu2ggFTree;
		TTree * VHTree;
		TTree * VHHadTree;
		TTree * VHLepTree;
		TTree * VBFTree;

		// Variables that will be stored
		Int_t tRun;
		Int_t tEvent;
		Int_t tlbn;

		Float_t tm4l_unconstrained;
		Float_t tm4lerr_unconstrained;
		Float_t tmZ1_unconstrained;
		Float_t tmZ2_unconstrained;

		Float_t tm4l_constrained;
		Float_t tm4lerr_constrained;
		Float_t tmZ1_constrained;
		Float_t tmZ2_constrained;

		Float_t tweight;
		Float_t tweight_corr;
		Float_t tweight_lumi;
		
		// Reco lepton information
		Float_t tZ1_lepplus_pt;
		Float_t tZ1_lepminus_pt;
		Float_t tZ2_lepplus_pt;
		Float_t tZ2_lepminus_pt;
		Float_t tZ1_lepplus_eta;
		Float_t tZ1_lepminus_eta;
		Float_t tZ2_lepplus_eta;
		Float_t tZ2_lepminus_eta;
		Float_t tZ1_lepplus_phi;
		Float_t tZ1_lepminus_phi;
		Float_t tZ2_lepplus_phi;
		Float_t tZ2_lepminus_phi;
		Float_t tZ1_lepplus_m;
		Float_t tZ1_lepminus_m;
		Float_t tZ2_lepplus_m;
		Float_t tZ2_lepminus_m;
		Float_t tZ1_lepplus_E;
		Float_t tZ1_lepminus_E;
		Float_t tZ2_lepplus_E;
		Float_t tZ2_lepminus_E;
		Int_t tZ1_lepplus_id;
		Int_t tZ1_lepminus_id;
		Int_t tZ2_lepplus_id;
		Int_t tZ2_lepminus_id;
		Float_t tZ1_lepplus_cov_mom;
		Float_t tZ1_lepminus_cov_mom;
		Float_t tZ2_lepplus_cov_mom;
		Float_t tZ2_lepminus_cov_mom;
		
		Float_t tZ1_lepplus_pt_uncorr_CB;
		Float_t tZ1_lepminus_pt_uncorr_CB;
		Float_t tZ2_lepplus_pt_uncorr_CB;
		Float_t tZ2_lepminus_pt_uncorr_CB;

		Float_t tZ1_lepplus_pt_uncorr_ID;
		Float_t tZ1_lepminus_pt_uncorr_ID;
		Float_t tZ2_lepplus_pt_uncorr_ID;
		Float_t tZ2_lepminus_pt_uncorr_ID;

		Float_t tZ1_lepplus_pt_uncorr_MS;
		Float_t tZ1_lepminus_pt_uncorr_MS;
		Float_t tZ2_lepplus_pt_uncorr_MS;
		Float_t tZ2_lepminus_pt_uncorr_MS;

		// true turth information
		Float_t tZ1_lepplus_pt_truth;
		Float_t tZ1_lepminus_pt_truth;
		Float_t tZ2_lepplus_pt_truth;
		Float_t tZ2_lepminus_pt_truth;
		Float_t tZ1_lepplus_eta_truth;
		Float_t tZ1_lepminus_eta_truth;
		Float_t tZ2_lepplus_eta_truth;
		Float_t tZ2_lepminus_eta_truth;
		Float_t tZ1_lepplus_phi_truth;
		Float_t tZ1_lepminus_phi_truth;
		Float_t tZ2_lepplus_phi_truth;
		Float_t tZ2_lepminus_phi_truth;
		Float_t tZ1_lepplus_m_truth;
		Float_t tZ1_lepminus_m_truth;
		Float_t tZ2_lepplus_m_truth;
		Float_t tZ2_lepminus_m_truth;
		
		// Truth information for reco matched
		Float_t tZ1_lepplus_pt_truth_bare;
		Float_t tZ1_lepminus_pt_truth_bare;
		Float_t tZ2_lepplus_pt_truth_bare;
		Float_t tZ2_lepminus_pt_truth_bare;
		Float_t tZ1_lepplus_eta_truth_bare;
		Float_t tZ1_lepminus_eta_truth_bare;
		Float_t tZ2_lepplus_eta_truth_bare;
		Float_t tZ2_lepminus_eta_truth_bare;
		Float_t tZ1_lepplus_phi_truth_bare;
		Float_t tZ1_lepminus_phi_truth_bare;
		Float_t tZ2_lepplus_phi_truth_bare;
		Float_t tZ2_lepminus_phi_truth_bare;
		Float_t tZ1_lepplus_m_truth_bare;
		Float_t tZ1_lepminus_m_truth_bare;
		Float_t tZ2_lepplus_m_truth_bare;
		Float_t tZ2_lepminus_m_truth_bare;
		
		// For LepSys
		Float_t Lep_cl_eta[4];
		Float_t Lep_E_ZeeStatUp[4];
		Float_t Lep_E_ZeeStatDown[4];
		Float_t Lep_E_ZeeSystUp[4];
		Float_t Lep_E_ZeeSystDown[4];
		Float_t Lep_E_ZeeAllUp[4];
		Float_t Lep_E_ZeeAllDown[4];
		Float_t Lep_E_PSUp[4];
		Float_t Lep_E_PSDown[4];
		Float_t Lep_E_S12Up[4];
		Float_t Lep_E_S12Down[4];
		Float_t Lep_E_MatIDUp[4];
		Float_t Lep_E_MatIDDown[4];
		Float_t Lep_E_MatCryoUp[4];
		Float_t Lep_E_MatCryoDown[4];
		Float_t Lep_E_MatCaloUp[4];
		Float_t Lep_E_MatCaloDown[4];
		Float_t Lep_E_LArCalibUp[4];
		Float_t Lep_E_LArCalibDown[4];
		Float_t Lep_E_LArUnconvCalibUp[4];
		Float_t Lep_E_LArUnconvCalibDown[4];
		Float_t Lep_E_LArElecUnconvUp[4];
		Float_t Lep_E_LArElecUnconvDown[4];
		Float_t Lep_E_LArElecCalibUp[4];
		Float_t Lep_E_LArElecCalibDown[4];
		Float_t Lep_E_GainUp[4];
		Float_t Lep_E_GainDown[4];
		Float_t Lep_E_G4Up[4];
		Float_t Lep_E_G4Down[4];
		Float_t Lep_E_MomentumUp[4];
		Float_t Lep_E_MomentumDown[4];
		Float_t Lep_E_ZSmearingUp[4];
		Float_t Lep_E_ZSmearingDown[4];
		Float_t Lep_E_SamplingTermUp[4];
		Float_t Lep_E_SamplingTermDown[4];
		Float_t Lep_E_MaterialIDUp[4];
		Float_t Lep_E_MaterialIDDown[4];
		Float_t Lep_E_MaterialCaloUp[4];
		Float_t Lep_E_MaterialCaloDown[4];
		Float_t Lep_E_MaterialGapUp[4];
		Float_t Lep_E_MaterialGapDown[4];
		Float_t Lep_E_MaterialCryoUp[4];
		Float_t Lep_E_MaterialCryoDown[4];
		Float_t Lep_E_PileUpUp[4];
		Float_t Lep_E_PileUpDown[4];
		Float_t Lep_E_Nom[4];

		Float_t Lep_E_MaterialAllUp[4];
		Float_t Lep_E_MaterialAllDown[4];
		

		//Split up
		Float_t Lep_E_PSUp_Barrel[4];
		Float_t Lep_E_PSDown_Barrel[4];
		Float_t Lep_E_S12Up_Barrel[4];
		Float_t Lep_E_S12Down_Barrel[4];
		Float_t Lep_E_MatCryoUp_Barrel[4];
		Float_t Lep_E_MatCryoDown_Barrel[4];
		Float_t Lep_E_MatCaloUp_Barrel[4];
		Float_t Lep_E_MatCaloDown_Barrel[4];
		Float_t Lep_E_LArUnconvCalibUp_Barrel[4];
		Float_t Lep_E_LArUnconvCalibDown_Barrel[4];
		Float_t Lep_E_LArElecUnconvUp_Barrel[4];
		Float_t Lep_E_LArElecUnconvDown_Barrel[4];
		Float_t Lep_E_LArCalibUp_Barrel[4];
		Float_t Lep_E_LArCalibDown_Barrel[4];
		
		Float_t Lep_E_PSUp_EC[4];
		Float_t Lep_E_PSDown_EC[4];
		Float_t Lep_E_S12Up_EC[4];
		Float_t Lep_E_S12Down_EC[4];
		Float_t Lep_E_MatCryoUp_EC[4];
		Float_t Lep_E_MatCryoDown_EC[4];
		Float_t Lep_E_MatCaloUp_EC[4];
		Float_t Lep_E_MatCaloDown_EC[4];
		Float_t Lep_E_LArUnconvCalibUp_EC[4];
		Float_t Lep_E_LArUnconvCalibDown_EC[4];
		Float_t Lep_E_LArElecUnconvUp_EC[4];
		Float_t Lep_E_LArElecUnconvDown_EC[4];
		Float_t Lep_E_LArCalibUp_EC[4];
		Float_t Lep_E_LArCalibDown_EC[4];

		Float_t Lep_E_MatIDUp1[4];
		Float_t Lep_E_MatIDDown1[4];
		Float_t Lep_E_MatIDUp2[4];
		Float_t Lep_E_MatIDDown2[4];
		Float_t Lep_E_MatIDUp3[4];
		Float_t Lep_E_MatIDDown3[4];
		Float_t Lep_E_MatIDUp4[4];
		Float_t Lep_E_MatIDDown4[4];




		// For weights
		Float_t tzVertWeight;
		Float_t tpileupWeight;
		Float_t tggFWeight;	
		Float_t tJHUWeight;				
		Float_t tttBarVeto;
		Float_t tzzOverlapVeto;
		Float_t tzBBVeto;
		Float_t ttrigEff;
		Float_t tlepEff;
		Float_t tmcEventWeight;
		Float_t tcrossSection;
		Float_t tbranchRatio;
		Float_t tlumi;


		Float_t tflagQuad;
};

#endif



