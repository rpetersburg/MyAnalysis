#include <stdlib.h>
#include <string>
#include "MyAnalysis/OutputTreeDi.h"


using namespace std;
////////////////////////////////////////////////////////////////////////////////////////
//				Constructor and Festructor
////////////////////////////////////////////////////////////////////////////////////////
OutputTreeDi::OutputTreeDi()
{
	// Clearing the vars
	clearVars();

	// Creating the trees
	mu2Tree = new TTree("tree_incl_2mu", "2mu");
	el2Tree = new TTree("tree_incl_2e", "2e");
	
	// Booking the trees
	bookTree(mu2Tree);
	bookTree(el2Tree);

}

OutputTreeDi::~OutputTreeDi()
{
	delete mu2Tree;
	delete el2Tree;
		
}

void OutputTreeDi::clearVars()
{
	tRun = -999;
	tEvent = -999;
	tlbn = -999;

	tm2l_unconstrained = -999;
	tm2lerr_unconstrained = -999;

	tm2l_constrained = -999;
	tm2lerr_constrained = -999;

	tweight = -999;
	tweight_corr = -999;
	tweight_lumi = -999;

	tZ1_lepplus_pt = -999;
	tZ1_lepminus_pt = -999;
	tZ1_lepplus_eta = -999;
	tZ1_lepminus_eta = -999;
	tZ1_lepplus_phi = -999;
	tZ1_lepminus_phi = -999;
	tZ1_lepplus_m = -999;
	tZ1_lepminus_m = -999;
	tZ1_lepplus_id = -999;
	tZ1_lepminus_id = -999;

	tZ1_lepplus_pt_truth = -999;
	tZ1_lepminus_pt_truth = -999;
	tZ1_lepplus_eta_truth = -999;
	tZ1_lepminus_eta_truth = -999;
	tZ1_lepplus_phi_truth = -999;
	tZ1_lepminus_phi_truth = -999;
	tZ1_lepplus_m_truth = -999;
	tZ1_lepminus_m_truth = -999;
	
	tZ1_lepplus_pt_truth_bare = -999;
	tZ1_lepminus_pt_truth_bare = -999;
	tZ1_lepplus_eta_truth_bare = -999;
	tZ1_lepminus_eta_truth_bare = -999;
	tZ1_lepplus_phi_truth_bare = -999;
	tZ1_lepminus_phi_truth_bare = -999;
	tZ1_lepplus_m_truth_bare = -999;
	tZ1_lepminus_m_truth_bare = -999;


	tZ1_lepplus_cov_mom = -999;
	tZ1_lepminus_cov_mom = -999;
	
	tzVertWeight = -999;
	tpileupWeight = -999;
	tggFWeight = -999;	
	tJHUWeight = -999;			
	tttBarVeto = -999;
	tzzOverlapVeto = -999;
	tzBBVeto = -999;
	ttrigEff = -999;
	tlepEff = -999;
	tmcEventWeight = -999;
	tcrossSection = -999;
	tbranchRatio = -999;
	tlumi = -999;

	//tflagQuad = -999;	
}

void OutputTreeDi::bookTree(TTree *tree)
{
	//cout<<"Booking Tree: "<<tree->getName() <<endl;

	tree->Branch("run" , 									&tRun , 								"run/I"					);		
	tree->Branch("event" , 									&tEvent ,                   			"event/I" 				);	
	tree->Branch("lbn" , 									&tlbn ,                     			"lbn/I" 				);	
	
	tree->Branch("m2l_unconstrained" , 						&tm2l_unconstrained ,       			"m2l_unconstrained/F"   );
	tree->Branch("m2lerr_unconstrained" , 					&tm2lerr_unconstrained ,    			"m2lerr_unconstrained/F"); 	
	tree->Branch("m2l_constrained" , 						&tm2l_constrained ,         			"m2l_constrained/F"     );
	tree->Branch("m2lerr_constrained" ,						&tm2lerr_constrained ,      			"m2lerr_constrained/F" 	);	
	
	tree->Branch("weight" , 								&tweight ,                  			"weight/F"  			);	
	tree->Branch("weight_corr" , 							&tweight_corr ,             			"weight_corr/F"  		);	
	tree->Branch("weight_lumi" , 							&tweight_lumi ,             			"weight_lumi/F"  		);	
	
	tree->Branch("Z1_lepplus_pt" , 							&tZ1_lepplus_pt ,           			"Z1_lepplus_pt/F"  		);	
	tree->Branch("Z1_lepminus_pt" , 						&tZ1_lepminus_pt ,          			"Z1_lepminus_pt/F"  	);	
	tree->Branch("Z1_lepplus_eta" , 						&tZ1_lepplus_eta ,          			"Z1_lepplus_eta/F"  	);	
	tree->Branch("Z1_lepminus_eta" , 						&tZ1_lepminus_eta ,         			"Z1_lepminus_eta/F"  	);	
	tree->Branch("Z1_lepplus_phi" , 						&tZ1_lepplus_phi ,          			"Z1_lepplus_phi/F"  	);	
	tree->Branch("Z1_lepminus_phi" ,						&tZ1_lepminus_phi ,         			"Z1_lepminus_phi/F" 	);	
	tree->Branch("Z1_lepplus_m" , 							&tZ1_lepplus_m ,            			"Z1_lepplus_m/F"  		);	
	tree->Branch("Z1_lepminus_m" , 							&tZ1_lepminus_m ,           			"Z1_lepminus_m/F"  		);	
	
	tree->Branch("Z1_lepplus_cov_mom" , 					&tZ1_lepplus_cov_mom ,      			"Z1_lepplus_cov_mom/F"   );	
	tree->Branch("Z1_lepminus_cov_mom" , 					&tZ1_lepminus_cov_mom ,     			"Z1_lepminus_cov_mom/F"  );	

	tree->Branch("Z1_lepplus_id" , 							&tZ1_lepplus_id ,            			"Z1_lepplus_id/I"  		);	
	tree->Branch("Z1_lepminus_id" , 						&tZ1_lepminus_id ,           			"Z1_lepminus_id/I"  	);	
	
	tree->Branch("Z1_lepplus_pt_uncorr_CB" , 				&tZ1_lepplus_pt_uncorr_CB ,           	"Z1_lepplus_pt_uncorr_CB/F"  	);	
	tree->Branch("Z1_lepminus_pt_uncorr_CB" , 				&tZ1_lepminus_pt_uncorr_CB ,          	"Z1_lepminus_pt_uncorr_CB/F"  	);	

	tree->Branch("Z1_lepplus_pt_uncorr_ID" , 				&tZ1_lepplus_pt_uncorr_ID ,           	"Z1_lepplus_pt_uncorr_ID/F"  	);	
	tree->Branch("Z1_lepminus_pt_uncorr_ID" , 				&tZ1_lepminus_pt_uncorr_ID ,          	"Z1_lepminus_pt_uncorr_ID/F"  	);	
	
	tree->Branch("Z1_lepplus_pt_uncorr_MS" , 				&tZ1_lepplus_pt_uncorr_MS ,           	"Z1_lepplus_pt_uncorr_MS/F"  	);	
	tree->Branch("Z1_lepminus_pt_uncorr_MS" , 				&tZ1_lepminus_pt_uncorr_MS ,          	"Z1_lepminus_pt_uncorr_MS/F"  	);	

	tree->Branch("Z1_lepplus_pt_truth_born" , 				&tZ1_lepplus_pt_truth ,     			"Z1_lepplus_pt_truth/F" 	);	
	tree->Branch("Z1_lepminus_pt_truth_born" , 				&tZ1_lepminus_pt_truth ,    			"Z1_lepminus_pt_truth/F"  	);	
	tree->Branch("Z1_lepplus_eta_truth_born" , 				&tZ1_lepplus_eta_truth ,    			"Z1_lepplus_eta_truth/F"  	);	
	tree->Branch("Z1_lepminus_eta_truth_born" , 			&tZ1_lepminus_eta_truth ,   			"Z1_lepminus_eta_truth/F"  	);	
	tree->Branch("Z1_lepplus_phi_truth_born" , 				&tZ1_lepplus_phi_truth ,    			"Z1_lepplus_phi_truth/F"  	);	
	tree->Branch("Z1_lepminus_phi_truth_born" ,				&tZ1_lepminus_phi_truth ,   			"Z1_lepminus_phi_truth/F" 	);	
	tree->Branch("Z1_lepplus_m_truth_born" , 				&tZ1_lepplus_m_truth ,      			"Z1_lepplus_m_truth/F"  	);	
	tree->Branch("Z1_lepminus_m_truth_born" , 				&tZ1_lepminus_m_truth ,     			"Z1_lepminus_m_truth/F"  	);	

	tree->Branch("Z1_lepplus_pt_truth_matched_bare",		&tZ1_lepplus_pt_truth_bare ,           	"Z1_lepplus_pt_truth_matched_bare/F"  	);	
	tree->Branch("Z1_lepminus_pt_truth_matched_bare" , 		&tZ1_lepminus_pt_truth_bare ,          	"Z1_lepminus_pt_truth_matched_bare/F"  	);	
	tree->Branch("Z1_lepplus_eta_truth_matched_bare" , 		&tZ1_lepplus_eta_truth_bare ,          	"Z1_lepplus_eta_truth_matched_bare/F"  	);	
	tree->Branch("Z1_lepminus_eta_truth_matched_bare" , 	&tZ1_lepminus_eta_truth_bare ,         	"Z1_lepminus_eta_truth_matched_bare/F" 	);	
	tree->Branch("Z1_lepplus_phi_truth_matched_bare" , 		&tZ1_lepplus_phi_truth_bare ,          	"Z1_lepplus_phi_truth_matched_bare/F"  	);	
	tree->Branch("Z1_lepminus_phi_truth_matched_bare" ,		&tZ1_lepminus_phi_truth_bare ,         	"Z1_lepminus_phi_truth_matched_bare/F" 	);	
	tree->Branch("Z1_lepplus_m_truth_matched_bare" , 		&tZ1_lepplus_m_truth_bare ,            	"Z1_lepplus_m_truth_matched_bare/F"  	);	
	tree->Branch("Z1_lepminus_m_truth_matched_bare" , 		&tZ1_lepminus_m_truth_bare ,           	"Z1_lepminus_m_truth_matched_bare/F"  	);	


	tree->Branch("zVertWeight" , 							&tzVertWeight ,             			"zVertWeight/F"  		);
	tree->Branch("ggFWeight" , 								&tggFWeight ,             				"ggFWeight/F"  			);	
	tree->Branch("JHUWeight" , 								&tJHUWeight ,             				"JHUWeight/F"  			);	
	tree->Branch("pileUpWeight" , 							&tpileupWeight ,            			"pileUpWeight/F"  		);	
	tree->Branch("ttBarVeto" , 								&tttBarVeto ,               			"ttBarVeto/F"  			);	
	tree->Branch("zzOverlapVeto" ,							&tzzOverlapVeto ,           			"zzOverlapVeto/F"  		);	
	tree->Branch("zBBVeto" , 								&tzBBVeto ,                 			"zBBVeto/F"  			);	
	tree->Branch("trigEff" , 								&ttrigEff ,            					"trigEff/F"  			);	
	tree->Branch("lepEff" , 								&tlepEff ,           					"lepEff/F"  			);
	tree->Branch("mcEventWeight" , 							&tmcEventWeight ,           			"mcEventWeight/F"  		);
	tree->Branch("crossSection" , 							&tcrossSection ,            			"crossSection/F"  		);	
	tree->Branch("branchRatio" , 							&tbranchRatio ,           				"branchRatio/F"  		);
	tree->Branch("lumi" , 									&tlumi ,           						"lumi/F"  				);

	//tree->Branch("flagQuad" , 								&tflagQuad ,           					"flagQuad/I"  			);	
}

void OutputTreeDi::fillTree(D3PDReader::Event *event, DiLepton * ZCan, Int_t type, Bool_t isMC)
{
	// Filling the variables
	if(isMC) tRun = event->eventinfo.mc_channel_number();
	else tRun = event->eventinfo.RunNumber();
	tEvent = event->eventinfo.EventNumber();;
	tlbn = event->eventinfo.lbn();

	tm2l_unconstrained 		= ZCan->mass/1000;
	tm2lerr_unconstrained 	= ZCan->massErr/1000;

	tm2l_constrained 		= ZCan->massZMassCons/1000;
	tm2lerr_constrained 	= ZCan->massErrZmassCons/1000;

	tweight 		= ZCan->weight;
	tweight_corr 	= ZCan->weight_corr;
	tweight_lumi 	= ZCan->weight_lumi;

	tZ1_lepplus_pt		= ZCan->getLepPlus()->get4Momentum()->Pt()/1000;
	tZ1_lepminus_pt 	= ZCan->getLepNeg()->get4Momentum()->Pt()/1000;

	tZ1_lepplus_eta 	= ZCan->getLepPlus()->get4Momentum()->Eta();	
	tZ1_lepminus_eta 	= ZCan->getLepNeg()->get4Momentum()->Eta();

	tZ1_lepplus_phi 	= ZCan->getLepPlus()->get4Momentum()->Phi();
	tZ1_lepminus_phi 	= ZCan->getLepNeg()->get4Momentum()->Phi();

	tZ1_lepplus_m 		= ZCan->getLepPlus()->lepMass/1000;
	tZ1_lepminus_m 		= ZCan->getLepNeg()->lepMass/1000;

	tZ1_lepplus_id 		= (Int_t) ZCan->getLepPlus()->lepID;
	tZ1_lepminus_id 	= (Int_t) ZCan->getLepNeg()->lepID;

	tZ1_lepplus_cov_mom  	= ZCan->getLepPlus()->covMomErr;
	tZ1_lepminus_cov_mom  	= ZCan->getLepNeg()->covMomErr;

	if(ZCan->getLepPlus()->m_momentumTruthRecoBare.Pt() && isMC)
	{
		tZ1_lepplus_pt_truth_bare 	= ZCan->getLepPlus()->m_momentumTruthRecoBare.Pt()/1000;
		tZ1_lepplus_eta_truth_bare 	= ZCan->getLepPlus()->m_momentumTruthRecoBare.Eta();
		tZ1_lepplus_phi_truth_bare 	= ZCan->getLepPlus()->m_momentumTruthRecoBare.Phi();
		tZ1_lepplus_m_truth_bare 	= ZCan->getLepPlus()->m_momentumTruthRecoBare.M()/1000;
	}
	if(ZCan->getLepNeg()->m_momentumTruthRecoBare.Pt() && isMC)
	{
		tZ1_lepminus_pt_truth_bare 	= ZCan->getLepNeg()->m_momentumTruthRecoBare.Pt()/1000;
		tZ1_lepminus_eta_truth_bare = ZCan->getLepNeg()->m_momentumTruthRecoBare.Eta();
		tZ1_lepminus_phi_truth_bare = ZCan->getLepNeg()->m_momentumTruthRecoBare.Phi();
		tZ1_lepminus_m_truth_bare 	= ZCan->getLepNeg()->m_momentumTruthRecoBare.M()/1000;
	}
	  
	if(ZCan->getLepPlus()->m_momentumTruthRecoBorn.Pt() && isMC)
	{
		tZ1_lepplus_pt_truth 	= ZCan->getLepPlus()->m_momentumTruthRecoBorn.Pt()/1000;
		tZ1_lepplus_eta_truth 	= ZCan->getLepPlus()->m_momentumTruthRecoBorn.Eta();
		tZ1_lepplus_phi_truth 	= ZCan->getLepPlus()->m_momentumTruthRecoBorn.Phi();
		tZ1_lepplus_m_truth 	= ZCan->getLepPlus()->m_momentumTruthRecoBorn.M()/1000;
	}
	if(ZCan->getLepNeg()->m_momentumTruthRecoBorn.Pt() && isMC)
	{
		tZ1_lepminus_pt_truth	= ZCan->getLepNeg()->m_momentumTruthRecoBorn.Pt()/1000;
		tZ1_lepminus_eta_truth 	= ZCan->getLepNeg()->m_momentumTruthRecoBorn.Eta();
		tZ1_lepminus_phi_truth 	= ZCan->getLepNeg()->m_momentumTruthRecoBorn.Phi();
		tZ1_lepminus_m_truth 	= ZCan->getLepNeg()->m_momentumTruthRecoBorn.M()/1000;
	}
	tZ1_lepplus_pt_uncorr_CB  	= ZCan->getLepPlus()->cb_pt_unsmeared/1000;
	tZ1_lepminus_pt_uncorr_CB  	= ZCan->getLepNeg()->cb_pt_unsmeared/1000;
	
	tZ1_lepplus_pt_uncorr_ID  	= ZCan->getLepPlus()->id_pt_unsmeared/1000;
	tZ1_lepminus_pt_uncorr_ID  	= ZCan->getLepNeg()->id_pt_unsmeared/1000;
	
	tZ1_lepplus_pt_uncorr_MS 	= ZCan->getLepPlus()->me_pt_unsmeared/1000;
	tZ1_lepminus_pt_uncorr_MS  	= ZCan->getLepNeg()->me_pt_unsmeared/1000;

	tzVertWeight 	= ZCan->zVertWeight;
	tpileupWeight 	= ZCan->pileupWeight;
	tggFWeight 		= ZCan->ggFWeight;
	tJHUWeight 		= ZCan->JHUWeight;	
	tttBarVeto 		= ZCan->ttBarVeto;
	tzzOverlapVeto 	= ZCan->zzOverlapVeto;
	tzBBVeto 		= ZCan->zBBVeto;
	ttrigEff 		= ZCan->trigEff;
	tlepEff 		= ZCan->lepEff;
	tmcEventWeight 	= ZCan->mcEventWeight;
	tcrossSection 	= ZCan->crossSection;
	tbranchRatio 	= ZCan->branchRatio;
	tlumi			= ZCan->lumi;
	//tflagQuad 		= higgs->flagQuad;



	// Filling the tree
	if(type == diLeptonType::_2mu) mu2Tree->Fill();
	else if(type == diLeptonType::_2e) el2Tree->Fill();
	else {cout<<"OutputTree: fillTree: Analysis type not recognized"<<endl;}


	// Clearing the vars for the next one
	clearVars();
}

void OutputTreeDi::saveTrees(TString filePath, TH1F* countingHist, TString sampleName)
{
	TFile *output = new TFile (filePath, "RECREATE");

	output->cd();
	// Counting Hist
	countingHist->Write();
	// Saving the trees
	mu2Tree->Write();
	el2Tree->Write();
	// Output Folder with the sample name for later checking
	output->mkdir(sampleName);
	cout<<"Sample Name: "<<sampleName<<endl;

	output->Close();


}


