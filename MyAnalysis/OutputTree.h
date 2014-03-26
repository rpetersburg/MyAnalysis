#ifndef OUTPUTTREE_H
#define OUTPUTTREE_H

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


class OutputTree 
{
	public :
		OutputTree();

		~OutputTree();
		
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

		Float_t tm4l_fsr;
		Float_t tm4lerr_fsr;		
		Float_t tmZ1_fsr;
		Float_t tmZ2_fsr;
		Int_t 	tfsrType;

		Float_t tm4l_constrained;
		Float_t tm4lerr_constrained;
		Float_t tmZ1_constrained;
		Float_t tmZ2_constrained;

		Float_t tm4l_unconstrained_ID;
		Float_t tm4lerr_unconstrained_ID;
		Float_t tm4l_constrained_ID;
		Float_t tm4lerr_constrained_ID;
		Float_t tm4l_unconstrained_MS;
		Float_t tm4lerr_unconstrained_MS;
		Float_t tm4l_constrained_MS;
		Float_t tm4lerr_constrained_MS;

		Float_t tweight;
		Float_t tweight_corr;
		Float_t tweight_lumi;
		
		Float_t tpt4l_unconstrained;
		Float_t ty4l_unconstrained;
		Float_t teta4l_unconstrained;
		Float_t tpt4l_fsr;
		Float_t ty4l_fsr;
		Float_t teta4l_fsr;	  
		Float_t tpt4l_constrained;
		Float_t ty4l_constrained;
		Float_t teta4l_constrained;	  
		
		Float_t tpt4l_truth_born;
		Float_t ty4l_truth_born;
		Float_t teta4l_truth_born;
		Float_t tphi4l_truth_born;
		Float_t tm4l_truth_born;

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
		Int_t tZ1_lepplus_id;
		Int_t tZ1_lepminus_id;
		Int_t tZ2_lepplus_id;
		Int_t tZ2_lepminus_id;

		Int_t tZ1_lepplus_turthParent;
		Int_t tZ1_lepminus_turthParent;
		Int_t tZ2_lepplus_turthParent;
		Int_t tZ2_lepminus_turthParent;

		Float_t tZ1_lepplus_cov_mom;
		Float_t tZ1_lepminus_cov_mom;
		Float_t tZ2_lepplus_cov_mom;
		Float_t tZ2_lepminus_cov_mom;
		
		Float_t tZ1_lepplus_pt_uncorr;
		Float_t tZ1_lepminus_pt_uncorr;
		Float_t tZ2_lepplus_pt_uncorr;
		Float_t tZ2_lepminus_pt_uncorr;

		Float_t tZ1_lepplus_pt_uncorr_ID;
		Float_t tZ1_lepminus_pt_uncorr_ID;
		Float_t tZ2_lepplus_pt_uncorr_ID;
		Float_t tZ2_lepminus_pt_uncorr_ID;

		Float_t tZ1_lepplus_pt_uncorr_MS;
		Float_t tZ1_lepminus_pt_uncorr_MS;
		Float_t tZ2_lepplus_pt_uncorr_MS;
		Float_t tZ2_lepminus_pt_uncorr_MS;

		Float_t tm4l_truth_matched_born;
		Float_t tmZ1_truth_matched_born;
		Float_t tmZ2_truth_matched_born;
		Float_t tm4l_truth_matched_bare;
		Float_t tmZ1_truth_matched_bare;
		Float_t tmZ2_truth_matched_bare;
		Float_t tm4l_truth_matched_dressed;
		Float_t tmZ1_truth_matched_dressed;
		Float_t tmZ2_truth_matched_dressed;

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

		Float_t tcthstr;
		Float_t tphi1;
		Float_t tcth1;
		Float_t tcth2;
		Float_t tphi;
		Float_t thiggspt;
		Float_t thiggseta;

		// For Jets
		Int_t	tn_jets;
		Float_t tdijet_invmass;
		Float_t tdijet_deltaeta;
		Float_t tleading_jet_pt;
		Float_t tleading_jet_eta;
		Float_t tleading_jet_phi;
		Float_t tleading_jet_m;
		Float_t tleading_jet_width;
		Float_t tleading_jet_nTrk;
		Float_t tsubleading_jet_pt;
		Float_t tsubleading_jet_eta;
		Float_t tsubleading_jet_phi;
		Float_t tsubleading_jet_m;
		Float_t tsubleading_jet_width;
		Float_t tsubleading_jet_nTrk;
		Float_t tthird_jet_pt;
		Float_t tthird_jet_eta;
		Float_t tthird_jet_phi;
		Float_t tthird_jet_m;
		Float_t tthird_jet_width;
		Float_t tthird_jet_nTrk;
		// Truth
		Int_t 	tn_jets_truth_bare;
		Float_t	tleading_jet_pt_truth_bare;

		// Fudicual Truth jets
		Int_t 	tn_jets_fid;
		Float_t	tleading_jet_pt_fid;
		Int_t 	tn_jets_truth_fid;
		Float_t	tleading_jet_pt_truth_fid;
		
		Float_t tBDT_discriminant_VBF;
		Float_t tBDT_discriminant_HadVH;
		
		Float_t tKD_discriminant;
		Float_t tBDT_discriminant;
		Float_t tBDTGuass_discriminant;
		Float_t tptSysupFac;
		Float_t tptSysdownFac;

		// For Categories
		Float_t tmissing_et;
		Float_t tleading_additional_lepton_pt;
		Float_t tleading_additional_lepton_eta;      
		Float_t tleading_additional_lepton_phi;     
		Float_t tleading_additional_lepton_m;
		Int_t tleading_additional_lepton_type;
		Int_t tleading_additional_lepton_type_truth_matched_bare;
		Float_t tsubleading_additional_lepton_pt;
		Float_t tsubleading_additional_lepton_eta;
		Float_t tsubleading_additional_lepton_phi;     
		Float_t tsubleading_additional_lepton_m;
		Int_t tsubleading_additional_lepton_type;
		Int_t tsubleading_additional_lepton_type_truth_matched_bare;

		Int_t tnpv;				
		Int_t tcalib;				
		Int_t teventType;

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

		// For ControlRegion
		Float_t tValTrackIso[4];
		Float_t tValCaloIso[4];
		Float_t tValD0Sign[4];
		Bool_t 	tTrackIso[4];
		Bool_t 	tCaloIso[4];
		Bool_t 	tD0Sign[4];
		Bool_t 	tElID[4];

		Float_t tflagQuad;

		Int_t tBCHCutMedium;
		Int_t tBCHCutTight;
		
};

#endif



