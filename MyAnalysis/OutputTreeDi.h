#ifndef OUTPUTTREEDI_H
#define OUTPUTTREEDI_H

#include <TTree.h>
#include <TString.h>
#include <TLorentzVector.h>

#include "D3PDReader/Event.h"

#include "MyAnalysis/DiLepton.h"
#include "MyAnalysis/macroDef.h"

#include <stdlib.h>
#include <math.h>
#include <iomanip>
#include <string>
#include <TH1F.h>


class OutputTreeDi 
{
	public :
		OutputTreeDi();

		~OutputTreeDi();
		
		void clearVars();
		
		void fillTree(D3PDReader::Event *event, DiLepton * ZCan, Int_t type, Bool_t isMC);
		void saveTrees(TString filePath, TH1F* countingHist, TString sampleName);
			
		void bookTree(TTree *tree);
	private :
		// TTree for Each channel
		TTree * mu2Tree;
		TTree * el2Tree;


		// Variables that will be stored
		Int_t tRun;
		Int_t tEvent;
		Int_t tlbn;

		Float_t tm2l_unconstrained;
		Float_t tm2lerr_unconstrained;
		
		Float_t tm2l_constrained;
		Float_t tm2lerr_constrained;

		Float_t tweight;
		Float_t tweight_corr;
		Float_t tweight_lumi;
		
		// Reco lepton information
		Float_t tZ1_lepplus_pt;
		Float_t tZ1_lepminus_pt;
		Float_t tZ1_lepplus_eta;
		Float_t tZ1_lepminus_eta;
		Float_t tZ1_lepplus_phi;
		Float_t tZ1_lepminus_phi;
		Float_t tZ1_lepplus_m;
		Float_t tZ1_lepminus_m;
		Int_t tZ1_lepplus_id;
		Int_t tZ1_lepminus_id;
		Float_t tZ1_lepplus_cov_mom;
		Float_t tZ1_lepminus_cov_mom;

		Float_t tZ1_lepplus_pt_uncorr_CB;
		Float_t tZ1_lepminus_pt_uncorr_CB;

		Float_t tZ1_lepplus_pt_uncorr_ID;
		Float_t tZ1_lepminus_pt_uncorr_ID;

		Float_t tZ1_lepplus_pt_uncorr_MS;
		Float_t tZ1_lepminus_pt_uncorr_MS;
		
		// true turth information
		Float_t tZ1_lepplus_pt_truth;
		Float_t tZ1_lepminus_pt_truth;
		Float_t tZ1_lepplus_eta_truth;
		Float_t tZ1_lepminus_eta_truth;
		Float_t tZ1_lepplus_phi_truth;
		Float_t tZ1_lepminus_phi_truth;
		Float_t tZ1_lepplus_m_truth;
		Float_t tZ1_lepminus_m_truth;
		
		// Truth information for reco matched
		Float_t tZ1_lepplus_pt_truth_bare;
		Float_t tZ1_lepminus_pt_truth_bare;
		Float_t tZ1_lepplus_eta_truth_bare;
		Float_t tZ1_lepminus_eta_truth_bare;
		Float_t tZ1_lepplus_phi_truth_bare;
		Float_t tZ1_lepminus_phi_truth_bare;
		Float_t tZ1_lepplus_m_truth_bare;
		Float_t tZ1_lepminus_m_truth_bare;
		
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

		//Float_t tflagQuad;
};

#endif



