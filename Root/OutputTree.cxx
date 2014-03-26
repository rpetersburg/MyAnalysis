#include <stdlib.h>
#include <string>
#include "MyAnalysis/OutputTree.h"


using namespace std;
////////////////////////////////////////////////////////////////////////////////////////
<<<<<<< HEAD
//				Constructor and Destructor
=======
//				Constructor and Festructor
>>>>>>> bc7b9ddaf72f0a41dfe1bb5d9068cc4b03444c0d
////////////////////////////////////////////////////////////////////////////////////////
OutputTree::OutputTree()
{
	// Clearing the vars
	clearVars();

	// Creating the trees
	mu4Tree = new TTree("tree_incl_4mu", "4mu");
	el4Tree = new TTree("tree_incl_4e", "4e");
	mu2el2Tree = new TTree("tree_incl_2mu2e", "2mu2e");
	el2mu2Tree = new TTree("tree_incl_2e2mu", "2e2mu");
	
	mu4ggFTree = new TTree("tree_ggF_4mu", "ggF:4mu");
	el4ggFTree = new TTree("tree_ggF_4e", "ggF:4e");
	mu2el2ggFTree = new TTree("tree_ggF_2mu2e", "ggF:2mu2e");
	el2mu2ggFTree = new TTree("tree_ggF_2e2mu", "ggF:2e2mu");
	
	VBFTree = new TTree("tree_VBF", "VBF:4l");
	
	VHTree = new TTree("tree_VH", "VH:4l");
	VHLepTree = new TTree("tree_VH_Lep", "VH:4l+lep");
	VHHadTree = new TTree("tree_VH_Had", "VH:4l+hadr");

	// Booking the trees
	bookTree(mu4Tree);
	bookTree(el4Tree);
	bookTree(mu2el2Tree);
	bookTree(el2mu2Tree);
	
	bookTree(mu4ggFTree);
	bookTree(el4ggFTree);
	bookTree(mu2el2ggFTree);
	bookTree(el2mu2ggFTree);
	
	bookTree(VBFTree);
	bookTree(VHTree);
	bookTree(VHLepTree);
	bookTree(VHHadTree);	
}

OutputTree::~OutputTree()
{
	delete mu4Tree;
	delete el4Tree;
	delete mu2el2Tree;
	delete el2mu2Tree;

	delete mu4ggFTree;
	delete el4ggFTree;
	delete mu2el2ggFTree;
	delete el2mu2ggFTree;

	delete VBFTree;
	delete VHTree;
	delete VHLepTree;
	delete VHHadTree;
	
}

void OutputTree::clearVars()
{
	tRun = -999;
	tEvent = -999;
	tlbn = -999;

	tm4l_unconstrained = -999;
	tm4lerr_unconstrained = -999;
	tmZ1_unconstrained = -999;
	tmZ2_unconstrained = -999;

	tm4l_fsr = -999;
	tm4lerr_fsr = -999;		
	tmZ1_fsr = -999;
	tmZ2_fsr = -999;
	tfsrType = -999;

	tm4l_constrained = -999;
	tm4lerr_constrained = -999;
	tmZ1_constrained = -999;
	tmZ2_constrained = -999;
	
	tm4l_unconstrained_ID = -999;
	tm4lerr_unconstrained_ID = -999;
	tm4l_constrained_ID = -999;
	tm4lerr_constrained_ID = -999;

	tm4l_unconstrained_MS = -999;
	tm4lerr_unconstrained_MS = -999;
	tm4l_constrained_MS = -999;
	tm4lerr_constrained_MS = -999;

	tweight = -999;
	tweight_corr = -999;
	tweight_lumi = -999;

	tm4l_truth_matched_born = -999;
	tmZ1_truth_matched_born = -999;
	tmZ2_truth_matched_born = -999;
	tm4l_truth_matched_bare = -999;
	tmZ1_truth_matched_bare = -999;
	tmZ2_truth_matched_bare = -999;
	tm4l_truth_matched_dressed = -999;
	tmZ1_truth_matched_dressed = -999;
	tmZ2_truth_matched_dressed = -999;

	tpt4l_unconstrained = -999;
	ty4l_unconstrained = -999;
	teta4l_unconstrained = -999;
	tpt4l_fsr = -999;
	ty4l_fsr = -999;
	teta4l_fsr = -999;	  
	tpt4l_constrained = -999;
	ty4l_constrained = -999;
	teta4l_constrained = -999;	  
	
	tpt4l_truth_born = -999;
	ty4l_truth_born = -999;
	teta4l_truth_born = -999;
<<<<<<< HEAD
=======
	tphi4l_truth_born = -999;
>>>>>>> bc7b9ddaf72f0a41dfe1bb5d9068cc4b03444c0d
	tm4l_truth_born = -999;

	tZ1_lepplus_pt = -999;
	tZ1_lepminus_pt = -999;
	tZ2_lepplus_pt = -999;
	tZ2_lepminus_pt = -999;
	tZ1_lepplus_eta = -999;
	tZ1_lepminus_eta = -999;
	tZ2_lepplus_eta = -999;
	tZ2_lepminus_eta = -999;
	tZ1_lepplus_phi = -999;
	tZ1_lepminus_phi = -999;
	tZ2_lepplus_phi = -999;
	tZ2_lepminus_phi = -999;
	tZ1_lepplus_m = -999;
	tZ1_lepminus_m = -999;
	tZ2_lepplus_m = -999;
	tZ2_lepminus_m = -999;
	tZ1_lepplus_id = -999;
	tZ1_lepminus_id = -999;
	tZ2_lepplus_id = -999;
	tZ2_lepminus_id = -999;
<<<<<<< HEAD
=======
	tZ1_lepplus_turthParent = -999;
	tZ1_lepminus_turthParent = -999;
	tZ2_lepplus_turthParent = -999;
	tZ2_lepminus_turthParent = -999;
>>>>>>> bc7b9ddaf72f0a41dfe1bb5d9068cc4b03444c0d

	tZ1_lepplus_pt_truth = -999;
	tZ1_lepminus_pt_truth = -999;
	tZ2_lepplus_pt_truth = -999;
	tZ2_lepminus_pt_truth = -999;
	tZ1_lepplus_eta_truth = -999;
	tZ1_lepminus_eta_truth = -999;
	tZ2_lepplus_eta_truth = -999;
	tZ2_lepminus_eta_truth = -999;
	tZ1_lepplus_phi_truth = -999;
	tZ1_lepminus_phi_truth = -999;
	tZ2_lepplus_phi_truth = -999;
	tZ2_lepminus_phi_truth = -999;
	tZ1_lepplus_m_truth = -999;
	tZ1_lepminus_m_truth = -999;
	tZ2_lepplus_m_truth = -999;
	tZ2_lepminus_m_truth = -999;
	
	tZ1_lepplus_pt_truth_bare = -999;
	tZ1_lepminus_pt_truth_bare = -999;
	tZ2_lepplus_pt_truth_bare = -999;
	tZ2_lepminus_pt_truth_bare = -999;
	tZ1_lepplus_eta_truth_bare = -999;
	tZ1_lepminus_eta_truth_bare = -999;
	tZ2_lepplus_eta_truth_bare = -999;
	tZ2_lepminus_eta_truth_bare = -999;
	tZ1_lepplus_phi_truth_bare = -999;
	tZ1_lepminus_phi_truth_bare = -999;
	tZ2_lepplus_phi_truth_bare = -999;
	tZ2_lepminus_phi_truth_bare = -999;
	tZ1_lepplus_m_truth_bare = -999;
	tZ1_lepminus_m_truth_bare = -999;
	tZ2_lepplus_m_truth_bare = -999;
	tZ2_lepminus_m_truth_bare = -999;


	tZ1_lepplus_cov_mom = -999;
	tZ1_lepminus_cov_mom = -999;
	tZ2_lepplus_cov_mom = -999;
	tZ2_lepminus_cov_mom = -999;
	
	tZ1_lepplus_pt_uncorr 	= -999;
	tZ1_lepminus_pt_uncorr 	= -999;
	tZ2_lepplus_pt_uncorr	= -999;
	tZ2_lepminus_pt_uncorr	= -999;

	tZ1_lepplus_pt_uncorr_ID = -999;
	tZ1_lepminus_pt_uncorr_ID = -999;
	tZ2_lepplus_pt_uncorr_ID = -999;
	tZ2_lepminus_pt_uncorr_ID = -999;

	tZ1_lepplus_pt_uncorr_MS = -999;
	tZ1_lepminus_pt_uncorr_MS = -999;
	tZ2_lepplus_pt_uncorr_MS = -999;
	tZ2_lepminus_pt_uncorr_MS = -999;

	tn_jets 				= -999;
	tdijet_invmass 			= -999;
	tdijet_deltaeta 		= -999;
	tleading_jet_pt 		= -999;
	tleading_jet_eta 		= -999;
	tleading_jet_phi		= -999;
	tleading_jet_m			= -999;
	tleading_jet_width		= -999;
	tleading_jet_nTrk		= -999;
	tsubleading_jet_pt		= -999;
	tsubleading_jet_eta		= -999;
	tsubleading_jet_phi		= -999;
	tsubleading_jet_m		= -999;
	tsubleading_jet_width	= -999;
	tsubleading_jet_nTrk	= -999;
	tthird_jet_pt			= -999;
	tthird_jet_eta			= -999;
	tthird_jet_phi			= -999;
	tthird_jet_m			= -999;
	tthird_jet_width		= -999;
	tthird_jet_nTrk			= -999;

	tBDT_discriminant_VBF = -999;
	tBDT_discriminant_HadVH = -999;
	tn_jets_truth_bare = -999;
	tleading_jet_pt_truth_bare = -999;
	
	tn_jets_fid = -999;
	tleading_jet_pt_fid = -999;
	tn_jets_truth_fid = -999;
	tleading_jet_pt_truth_fid = -999;

	tKD_discriminant = -999;
	tBDT_discriminant = -999;
	tBDTGuass_discriminant = -999;
	tptSysupFac = -999;
	tptSysdownFac = -999;

	tmissing_et = -999;
	tleading_additional_lepton_pt 		= -999;
	tleading_additional_lepton_eta 		= -999;
	tleading_additional_lepton_phi 		= -999;
	tleading_additional_lepton_m 		= -999;
	tleading_additional_lepton_type		= -999;
	tleading_additional_lepton_type_truth_matched_bare = -999;
	tsubleading_additional_lepton_pt 	= -999;
	tsubleading_additional_lepton_eta 	= -999;
	tsubleading_additional_lepton_phi 	= -999;
	tsubleading_additional_lepton_m		= -999;
	tsubleading_additional_lepton_type	= -999;
	tsubleading_additional_lepton_type_truth_matched_bare = -999;

	tnpv = -999;
	tcalib = -999;
	teventType = -999;

	tcthstr = -999;
	tphi1 = -999;
	tcth1 = -999;
	tcth2 = -999;
	tphi = -999;
	thiggspt = -999;
	thiggseta = -999;

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
	// for CR
	for(Int_t i = 0; i < 4; i++)
	{
		tValTrackIso[i] = -999; 
		tValCaloIso[i] 	= -999; 
		tValD0Sign[i] 	= -999; 

		tTrackIso[i]	= -999; 
		tCaloIso[i] 	= -999; 
		tD0Sign[i] 		= -999;
		tElID[i] 			= -999;
	}
	tflagQuad = -999;	
<<<<<<< HEAD
}

=======

	tBCHCutMedium = -999;
	tBCHCutTight = -999;
}



>>>>>>> bc7b9ddaf72f0a41dfe1bb5d9068cc4b03444c0d
void OutputTree::bookTree(TTree *tree)
{
	//cout<<"Booking Tree: "<<tree->getName() <<endl;

	tree->Branch("run" , 									&tRun , 								"run/I"					);		
	tree->Branch("event" , 									&tEvent ,                   			"event/I" 				);	
	tree->Branch("lbn" , 									&tlbn ,                     			"lbn/I" 				);	
	
	tree->Branch("m4l_unconstrained" , 						&tm4l_unconstrained ,       			"m4l_unconstrained/F"   );
	tree->Branch("m4lerr_unconstrained" , 					&tm4lerr_unconstrained ,    			"m4lerr_unconstrained/F"); 	
	tree->Branch("mZ1_unconstrained" , 						&tmZ1_unconstrained ,       			"mZ1_unconstrained/F"	);	
	tree->Branch("mZ2_unconstrained" , 						&tmZ2_unconstrained ,       			"mZ2_unconstrained/F" 	);	
	
	tree->Branch("m4l_fsr" , 								&tm4l_fsr ,                 			"m4l_fsr/F" 			);
	tree->Branch("mZ1_fsr" ,			 					&tmZ1_fsr ,                 			"mZ1_fsr/F"			 	);
	tree->Branch("mZ2_fsr" ,			 					&tmZ2_fsr ,                 			"mZ2_fsr/F"			 	);
	tree->Branch("m4lerr_fsr" ,			 					&tm4lerr_fsr ,              			"m4lerr_fsr/F"			);	
	tree->Branch("fsr_type" ,			 					&tfsrType ,              				"fsr_type/I"			);	
	
	tree->Branch("m4l_constrained" , 						&tm4l_constrained ,         			"m4l_constrained/F"     );
	tree->Branch("m4lerr_constrained" ,						&tm4lerr_constrained ,      			"m4lerr_constrained/F" 	);	
	tree->Branch("mZ1_constrained" , 						&tmZ1_constrained ,         			"mZ1_constrained/F"  	);	
	tree->Branch("mZ2_constrained" , 						&tmZ2_constrained ,         			"mZ2_constrained/F"  	);	

	tree->Branch("m4l_unconstrained_MS" , 					&tm4l_unconstrained_MS ,       			"m4l_unconstrained_MS/F"   );
	tree->Branch("m4lerr_unconstrained_MS" , 				&tm4lerr_unconstrained_MS ,    			"m4lerr_unconstrained_MS/F");
	tree->Branch("m4l_constrained_MS" , 					&tm4l_constrained_MS ,         			"m4l_constrained_MS/F"     );
	tree->Branch("m4lerr_constrained_MS" ,					&tm4lerr_constrained_MS ,      			"m4lerr_constrained_MS/F"  );
	tree->Branch("m4l_unconstrained_ID" , 					&tm4l_unconstrained_ID ,       			"m4l_unconstrained_ID/F"   );
	tree->Branch("m4lerr_unconstrained_ID" , 				&tm4lerr_unconstrained_ID ,    			"m4lerr_unconstrained_ID/F");
	tree->Branch("m4l_constrained_ID" , 					&tm4l_constrained_ID ,         			"m4l_constrained_ID/F"     );
	tree->Branch("m4lerr_constrained_ID" ,					&tm4lerr_constrained_ID ,      			"m4lerr_constrained_ID/F"  );

	tree->Branch("weight" , 								&tweight ,                  			"weight/F"  			);	
	tree->Branch("weight_corr" , 							&tweight_corr ,             			"weight_corr/F"  		);	
	tree->Branch("weight_lumi" , 							&tweight_lumi ,             			"weight_lumi/F"  		);	
	
	tree->Branch("pt4l_unconstrained" , 					&tpt4l_unconstrained ,       			"pt4l_unconstrained/F"  	); 
	tree->Branch("y4l_unconstrained" , 						&ty4l_unconstrained ,       			"y4l_unconstrained/F"  	); 
	tree->Branch("eta4l_unconstrained" , 					&teta4l_unconstrained ,       			"eta4l_unconstrained/F"  	); 
	tree->Branch("pt4l_fsr" , 								&tpt4l_fsr ,       						"pt4l_fsr/F"  	); 
	tree->Branch("y4l_fsr" , 								&ty4l_fsr ,       						"y4l_fsr/F"  	); 
	tree->Branch("eta4l_fsr" , 								&teta4l_fsr ,       					"eta4l_fsr/F"  	); 
	tree->Branch("pt4l_constrained" , 						&tpt4l_constrained ,       				"pt4l_constrained/F"  	); 
	tree->Branch("y4l_constrained" , 						&ty4l_constrained ,       				"y4l_constrained/F"  	); 
	tree->Branch("eta4l_constrained" , 						&teta4l_constrained ,       			"eta4l_constrained/F"  	); 
	
	tree->Branch("pt4l_truth_born" , 						&tpt4l_truth_born ,       				"pt4l_truth_born/F"  	); 
	tree->Branch("y4l_truth_born" , 						&ty4l_truth_born ,       				"y4l_truth_born/F"  	); 
	tree->Branch("eta4l_truth_born" , 						&teta4l_truth_born ,       				"eta4l_truth_born/F"  	); 
<<<<<<< HEAD
=======
	tree->Branch("phi4l_truth_born" , 						&tphi4l_truth_born ,       				"phi4l_truth_born/F"  	); 
>>>>>>> bc7b9ddaf72f0a41dfe1bb5d9068cc4b03444c0d
	tree->Branch("m4l_truth_born" , 						&tm4l_truth_born ,       				"m4l_truth_born/F"  	);

	tree->Branch("Z1_lepplus_pt" , 							&tZ1_lepplus_pt ,           			"Z1_lepplus_pt/F"  		);	
	tree->Branch("Z1_lepminus_pt" , 						&tZ1_lepminus_pt ,          			"Z1_lepminus_pt/F"  	);	
	tree->Branch("Z2_lepplus_pt" , 							&tZ2_lepplus_pt ,           			"Z2_lepplus_pt/F"  		);	
	tree->Branch("Z2_lepminus_pt" , 						&tZ2_lepminus_pt ,          			"Z2_lepminus_pt/F"  	);	
	tree->Branch("Z1_lepplus_eta" , 						&tZ1_lepplus_eta ,          			"Z1_lepplus_eta/F"  	);	
	tree->Branch("Z1_lepminus_eta" , 						&tZ1_lepminus_eta ,         			"Z1_lepminus_eta/F"  	);	
	tree->Branch("Z2_lepplus_eta" , 						&tZ2_lepplus_eta ,          			"Z2_lepplus_eta/F"  	);	
	tree->Branch("Z2_lepminus_eta" , 						&tZ2_lepminus_eta ,         			"Z2_lepminus_eta/F"  	);	
	tree->Branch("Z1_lepplus_phi" , 						&tZ1_lepplus_phi ,          			"Z1_lepplus_phi/F"  	);	
	tree->Branch("Z1_lepminus_phi" ,						&tZ1_lepminus_phi ,         			"Z1_lepminus_phi/F" 	);	
	tree->Branch("Z2_lepplus_phi" , 						&tZ2_lepplus_phi ,          			"Z2_lepplus_phi/F"  	);	
	tree->Branch("Z2_lepminus_phi" , 						&tZ2_lepminus_phi ,         			"Z2_lepminus_phi/F"  	);	
	tree->Branch("Z1_lepplus_m" , 							&tZ1_lepplus_m ,            			"Z1_lepplus_m/F"  		);	
	tree->Branch("Z1_lepminus_m" , 							&tZ1_lepminus_m ,           			"Z1_lepminus_m/F"  		);	
	tree->Branch("Z2_lepplus_m" , 							&tZ2_lepplus_m ,            			"Z2_lepplus_m/F"  		);	
	tree->Branch("Z2_lepminus_m" , 							&tZ2_lepminus_m ,           			"Z2_lepminus_m/F"  		);
	
	tree->Branch("Z1_lepplus_cov_mom" , 					&tZ1_lepplus_cov_mom ,      			"Z1_lepplus_cov_mom/F"   );	
	tree->Branch("Z1_lepminus_cov_mom" , 					&tZ1_lepminus_cov_mom ,     			"Z1_lepminus_cov_mom/F"  );	
	tree->Branch("Z2_lepplus_cov_mom" , 					&tZ2_lepplus_cov_mom ,      			"Z2_lepplus_cov_mom/F"   );	
	tree->Branch("Z2_lepminus_cov_mom" , 					&tZ2_lepminus_cov_mom ,     			"Z2_lepminus_cov_mom/F"  );	

	tree->Branch("Z1_lepplus_id" , 							&tZ1_lepplus_id ,            			"Z1_lepplus_id/I"  		);	
	tree->Branch("Z1_lepminus_id" , 						&tZ1_lepminus_id ,           			"Z1_lepminus_id/I"  	);	
	tree->Branch("Z2_lepplus_id" , 							&tZ2_lepplus_id ,            			"Z2_lepplus_id/I"  		);	
	tree->Branch("Z2_lepminus_id" , 						&tZ2_lepminus_id ,           			"Z2_lepminus_id/I"  	);

<<<<<<< HEAD
=======
	tree->Branch("Z1_lepplus_turthParent" , 				&tZ1_lepplus_turthParent ,           	"Z1_lepplus_turthParent/I"  	);	
	tree->Branch("Z1_lepminus_turthParent" , 				&tZ1_lepminus_turthParent ,    			"Z1_lepminus_turthParent/I"  	);	
	tree->Branch("Z2_lepplus_turthParent" , 				&tZ2_lepplus_turthParent ,     			"Z2_lepplus_turthParent/I"  	);	
	tree->Branch("Z2_lepminus_turthParent" , 				&tZ2_lepminus_turthParent ,    			"Z2_lepminus_turthParent/I"  	);

>>>>>>> bc7b9ddaf72f0a41dfe1bb5d9068cc4b03444c0d
	tree->Branch("Z1_lepplus_pt_uncorr" , 					&tZ1_lepplus_pt_uncorr ,  	         	"Z1_lepplus_pt_uncorr/F"  	);	
	tree->Branch("Z1_lepminus_pt_uncorr" , 					&tZ1_lepminus_pt_uncorr , 	         	"Z1_lepminus_pt_uncorr/F"  	);	
	tree->Branch("Z2_lepplus_pt_uncorr" , 					&tZ2_lepplus_pt_uncorr ,  	         	"Z2_lepplus_pt_uncorr/F" 	);	
	tree->Branch("Z2_lepminus_pt_uncorr" , 					&tZ2_lepminus_pt_uncorr , 	         	"Z2_lepminus_pt_uncorr/F"  	);

	tree->Branch("Z1_lepplus_pt_uncorr_ID" , 				&tZ1_lepplus_pt_uncorr_ID ,           	"Z1_lepplus_pt_uncorr_ID/F"  	);	
	tree->Branch("Z1_lepminus_pt_uncorr_ID" , 				&tZ1_lepminus_pt_uncorr_ID ,          	"Z1_lepminus_pt_uncorr_ID/F"  	);	
	tree->Branch("Z2_lepplus_pt_uncorr_ID" , 				&tZ2_lepplus_pt_uncorr_ID ,           	"Z2_lepplus_pt_uncorr_ID/F" 	);	
	tree->Branch("Z2_lepminus_pt_uncorr_ID" , 				&tZ2_lepminus_pt_uncorr_ID ,          	"Z2_lepminus_pt_uncorr_ID/F"  	);
	
	tree->Branch("Z1_lepplus_pt_uncorr_MS" , 				&tZ1_lepplus_pt_uncorr_MS ,           	"Z1_lepplus_pt_uncorr_MS/F"  	);	
	tree->Branch("Z1_lepminus_pt_uncorr_MS" , 				&tZ1_lepminus_pt_uncorr_MS ,          	"Z1_lepminus_pt_uncorr_MS/F"  	);	
	tree->Branch("Z2_lepplus_pt_uncorr_MS" , 				&tZ2_lepplus_pt_uncorr_MS ,           	"Z2_lepplus_pt_uncorr_MS/F" 	);	
	tree->Branch("Z2_lepminus_pt_uncorr_MS" , 				&tZ2_lepminus_pt_uncorr_MS ,          	"Z2_lepminus_pt_uncorr_MS/F"  	);

	tree->Branch("m4l_truth_matched_born" , 				&tm4l_truth_matched_born ,          	"m4l_truth_matched_born/F"  	);	
	tree->Branch("mZ1_truth_matched_born" , 				&tmZ1_truth_matched_born ,          	"mZ1_truth_matched_born/F"  	);	
	tree->Branch("mZ2_truth_matched_born" , 				&tmZ2_truth_matched_born ,          	"mZ2_truth_matched_born/F"  	);	
	tree->Branch("m4l_truth_matched_bare" , 				&tm4l_truth_matched_bare ,          	"m4l_truth_matched_bare/F"  	);	
	tree->Branch("mZ1_truth_matched_bare" , 				&tmZ1_truth_matched_bare ,          	"mZ1_truth_matched_bare/F"  	);	
	tree->Branch("mZ2_truth_matched_bare" , 				&tmZ2_truth_matched_bare ,          	"mZ2_truth_matched_bare/F"  	);	
    tree->Branch("m4l_truth_matched_dressed" , 				&tm4l_truth_matched_dressed ,       	"m4l_truth_matched_dressed/F"  	);	
	tree->Branch("mZ1_truth_matched_dressed" , 				&tmZ1_truth_matched_dressed ,       	"mZ1_truth_matched_dressed/F"  	);	
	tree->Branch("mZ2_truth_matched_dressed" , 				&tmZ2_truth_matched_dressed ,       	"mZ2_truth_matched_dressed/F"  	);

	tree->Branch("Z1_lepplus_pt_truth_born" , 				&tZ1_lepplus_pt_truth ,     			"Z1_lepplus_pt_truth/F" 	);	
	tree->Branch("Z1_lepminus_pt_truth_born" , 				&tZ1_lepminus_pt_truth ,    			"Z1_lepminus_pt_truth/F"  	);	
	tree->Branch("Z2_lepplus_pt_truth_born" , 				&tZ2_lepplus_pt_truth ,     			"Z2_lepplus_pt_truth/F"  	);	
	tree->Branch("Z2_lepminus_pt_truth_born" , 				&tZ2_lepminus_pt_truth ,    			"Z2_lepminus_pt_truth/F"  	);	
	tree->Branch("Z1_lepplus_eta_truth_born" , 				&tZ1_lepplus_eta_truth ,    			"Z1_lepplus_eta_truth/F"  	);	
	tree->Branch("Z1_lepminus_eta_truth_born" , 			&tZ1_lepminus_eta_truth ,   			"Z1_lepminus_eta_truth/F"  	);	
	tree->Branch("Z2_lepplus_eta_truth_born" , 				&tZ2_lepplus_eta_truth ,    			"Z2_lepplus_eta_truth/F"  	);	
	tree->Branch("Z2_lepminus_eta_truth_born" , 			&tZ2_lepminus_eta_truth ,   			"Z2_lepminus_eta_truth/F"  	);	
	tree->Branch("Z1_lepplus_phi_truth_born" , 				&tZ1_lepplus_phi_truth ,    			"Z1_lepplus_phi_truth/F"  	);	
	tree->Branch("Z1_lepminus_phi_truth_born" ,				&tZ1_lepminus_phi_truth ,   			"Z1_lepminus_phi_truth/F" 	);	
	tree->Branch("Z2_lepplus_phi_truth_born" , 				&tZ2_lepplus_phi_truth ,    			"Z2_lepplus_phi_truth/F"  	);	
	tree->Branch("Z2_lepminus_phi_truth_born" , 			&tZ2_lepminus_phi_truth ,   			"Z2_lepminus_phi_truth/F"  	);	
	tree->Branch("Z1_lepplus_m_truth_born" , 				&tZ1_lepplus_m_truth ,      			"Z1_lepplus_m_truth/F"  	);	
	tree->Branch("Z1_lepminus_m_truth_born" , 				&tZ1_lepminus_m_truth ,     			"Z1_lepminus_m_truth/F"  	);	
	tree->Branch("Z2_lepplus_m_truth_born" , 				&tZ2_lepplus_m_truth ,      			"Z2_lepplus_m_truth/F"  	);	
	tree->Branch("Z2_lepminus_m_truth_born" , 				&tZ2_lepminus_m_truth ,     			"Z2_lepminus_m_truth/F"  	);

	tree->Branch("Z1_lepplus_pt_truth_matched_bare",		&tZ1_lepplus_pt_truth_bare ,           	"Z1_lepplus_pt_truth_matched_bare/F"  	);	
	tree->Branch("Z1_lepminus_pt_truth_matched_bare" , 		&tZ1_lepminus_pt_truth_bare ,          	"Z1_lepminus_pt_truth_matched_bare/F"  	);	
	tree->Branch("Z2_lepplus_pt_truth_matched_bare" , 		&tZ2_lepplus_pt_truth_bare ,           	"Z2_lepplus_pt_truth_matched_bare/F"  	);	
	tree->Branch("Z2_lepminus_pt_truth_matched_bare" , 		&tZ2_lepminus_pt_truth_bare ,          	"Z2_lepminus_pt_truth_matched_bare/F"  	);	
	tree->Branch("Z1_lepplus_eta_truth_matched_bare" , 		&tZ1_lepplus_eta_truth_bare ,          	"Z1_lepplus_eta_truth_matched_bare/F"  	);	
	tree->Branch("Z1_lepminus_eta_truth_matched_bare" , 	&tZ1_lepminus_eta_truth_bare ,         	"Z1_lepminus_eta_truth_matched_bare/F" 	);	
	tree->Branch("Z2_lepplus_eta_truth_matched_bare" , 		&tZ2_lepplus_eta_truth_bare ,          	"Z2_lepplus_eta_truth_matched_bare/F"  	);	
	tree->Branch("Z2_lepminus_eta_truth_matched_bare" , 	&tZ2_lepminus_eta_truth_bare ,         	"Z2_lepminus_eta_truth_matched_bare/F" 	);	
	tree->Branch("Z1_lepplus_phi_truth_matched_bare" , 		&tZ1_lepplus_phi_truth_bare ,          	"Z1_lepplus_phi_truth_matched_bare/F"  	);	
	tree->Branch("Z1_lepminus_phi_truth_matched_bare" ,		&tZ1_lepminus_phi_truth_bare ,         	"Z1_lepminus_phi_truth_matched_bare/F" 	);	
	tree->Branch("Z2_lepplus_phi_truth_matched_bare" , 		&tZ2_lepplus_phi_truth_bare ,          	"Z2_lepplus_phi_truth_matched_bare/F"  	);	
	tree->Branch("Z2_lepminus_phi_truth_matched_bare" , 	&tZ2_lepminus_phi_truth_bare ,         	"Z2_lepminus_phi_truth_matched_bare/F" 	);	
	tree->Branch("Z1_lepplus_m_truth_matched_bare" , 		&tZ1_lepplus_m_truth_bare ,            	"Z1_lepplus_m_truth_matched_bare/F"  	);	
	tree->Branch("Z1_lepminus_m_truth_matched_bare" , 		&tZ1_lepminus_m_truth_bare ,           	"Z1_lepminus_m_truth_matched_bare/F"  	);	
	tree->Branch("Z2_lepplus_m_truth_matched_bare" , 		&tZ2_lepplus_m_truth_bare ,            	"Z2_lepplus_m_truth_matched_bare/F"  	);	
	tree->Branch("Z2_lepminus_m_truth_matched_bare" , 		&tZ2_lepminus_m_truth_bare ,           	"Z2_lepminus_m_truth_matched_bare/F"  	);	

	tree->Branch("cthstr" , 								&tcthstr ,                  			"cthstr/F"  			);	
	tree->Branch("phi1" , 									&tphi1 ,                    			"phi1/F"  				);	
	tree->Branch("cth1" , 									&tcth1 ,                    			"cth1/F"  				);	
	tree->Branch("cth2" , 									&tcth2 ,                    			"cth2/F"  				);	
	tree->Branch("phi" , 									&tphi ,                     			"phi/F"  				);	
	
	tree->Branch("npv" , 									&tnpv ,      							"npv/I"   );		
	tree->Branch("data_calib_type" , 						&tcalib ,      							"data_calib_type/I"   );		
	tree->Branch("event_type" , 							&teventType ,      						"event_type/I"   );		

	tree->Branch("n_jets" , 								&tn_jets ,      						"n_jets/I"   );		
	tree->Branch("dijet_invmass" , 							&tdijet_invmass ,      					"dijet_invmass/F"   );	
	tree->Branch("dijet_deltaeta" , 						&tdijet_deltaeta ,     					"dijet_deltaeta/F"  );	
	
	tree->Branch("leading_jet_pt" , 						&tleading_jet_pt ,      				"leading_jet_pt/F"   );	
	tree->Branch("leading_jet_eta" , 						&tleading_jet_eta ,     				"leading_jet_eta/F"  );
	tree->Branch("leading_jet_phi" , 						&tleading_jet_phi ,     				"leading_jet_phi/F"  );
	tree->Branch("leading_jet_m" , 							&tleading_jet_m ,	     				"leading_jet_m/F"  	 );
	tree->Branch("leading_jet_width" , 						&tleading_jet_width ,	     			"leading_jet_width/F");
	tree->Branch("leading_jet_nTrk" , 						&tleading_jet_nTrk ,	     			"leading_jet_nTrk/F" );
	
	tree->Branch("subleading_jet_pt" , 						&tsubleading_jet_pt ,      				"subleading_jet_pt/F"   );	
	tree->Branch("subleading_jet_eta" , 					&tsubleading_jet_eta ,     				"subleading_jet_eta/F"  );
	tree->Branch("subleading_jet_phi" , 					&tsubleading_jet_phi ,     				"subleading_jet_phi/F"  );
	tree->Branch("subleading_jet_m" , 						&tsubleading_jet_m ,	   				"subleading_jet_m/F"  	 );
	tree->Branch("subleading_jet_width" , 					&tsubleading_jet_width ,	     		"subleading_jet_width/F");
	tree->Branch("subleading_jet_nTrk" , 					&tsubleading_jet_nTrk ,	     			"subleading_jet_nTrk/F" );
	
	tree->Branch("third_jet_pt" , 							&tthird_jet_pt ,      					"third_jet_pt/F"   );	
	tree->Branch("third_jet_eta" , 							&tthird_jet_eta ,     					"third_jet_eta/F"  );
	tree->Branch("third_jet_phi" , 							&tthird_jet_phi ,     					"third_jet_phi/F"  );
	tree->Branch("third_jet_m" , 							&tthird_jet_m ,	     					"third_jet_m/F"  	 );
	tree->Branch("third_jet_width" , 						&tthird_jet_width ,	     				"third_jet_width/F");
	tree->Branch("third_jet_nTrk" , 						&tthird_jet_nTrk ,	     				"third_jet_nTrk/F" );

	tree->Branch("BDT_discriminant_VBF" , 					&tBDT_discriminant_VBF ,     			"BDT_discriminant_VBF/F"   );	
	tree->Branch("BDT_discriminant_HadVH" , 				&tBDT_discriminant_HadVH ,   			"BDT_discriminant_HadVH/F"  );
	tree->Branch("n_jets_truth_bare" , 						&tn_jets_truth_bare ,      				"n_jets_truth_bare/I"   );	
	tree->Branch("leading_jet_pt_truth_bare" , 				&tleading_jet_pt_truth_bare ,      		"leading_jet_pt_truth_bare/F");	

	tree->Branch("KD_discriminant" , 						&tKD_discriminant ,     				"KD_discriminant/F"   );	
	tree->Branch("BDT_discriminant" , 						&tBDT_discriminant ,   					"BDT_discriminant/F"  );
	tree->Branch("outBDT_gauss" , 							&tBDTGuass_discriminant ,     			"outBDT_gauss/F"   );	
	tree->Branch("ptSysupFac" , 							&tptSysupFac ,   						"ptSysupFac/F"  );
	tree->Branch("ptSysdownFac" , 							&tptSysdownFac ,   						"ptSysdownFac/F"  );

	tree->Branch("missing_et" , 							&tmissing_et ,   						"missing_et/F"  );
	tree->Branch("leading_additional_lepton_pt" , 			&tleading_additional_lepton_pt ,   		"leading_additional_lepton_pt/F" 	 	);
	tree->Branch("leading_additional_lepton_eta" , 			&tleading_additional_lepton_eta ,   	"leading_additional_lepton_eta/F" 	 	);
	tree->Branch("leading_additional_lepton_phi" , 			&tleading_additional_lepton_phi ,   	"leading_additional_lepton_phi/F"		);
	tree->Branch("leading_additional_lepton_m" , 			&tleading_additional_lepton_m ,   		"leading_additional_lepton_m/F"  		);
	tree->Branch("leading_additional_lepton_type" , 		&tleading_additional_lepton_type ,   	"leading_additional_lepton_type/I"  	);
	tree->Branch("leading_additional_lepton_type_truth_matched_bare" , &tleading_additional_lepton_type_truth_matched_bare ,"leading_additional_lepton_type_truth_matched_bare/I");
	tree->Branch("subleading_additional_lepton_pt" , 		&tsubleading_additional_lepton_pt ,   	"subleading_additional_lepton_pt/F"  	);
	tree->Branch("subleading_additional_lepton_eta" , 		&tsubleading_additional_lepton_eta ,   	"subleading_additional_lepton_eta/F"  	);
	tree->Branch("subleading_additional_lepton_phi" , 		&tsubleading_additional_lepton_phi ,   	"subleading_additional_lepton_phi/F"  	);
	tree->Branch("subleading_additional_lepton_m" , 		&tsubleading_additional_lepton_m ,   	"subleading_additional_lepton_m/F"  	);
	tree->Branch("subleading_additional_lepton_type" , 		&tsubleading_additional_lepton_type ,  	"subleading_additional_lepton_type/I"  	);
	tree->Branch("subleading_additional_lepton_type_truth_matched_bare" , &tsubleading_additional_lepton_type_truth_matched_bare ,"subleading_additional_lepton_type_truth_matched_bare/I");

<<<<<<< HEAD

=======
>>>>>>> bc7b9ddaf72f0a41dfe1bb5d9068cc4b03444c0d
	tree->Branch("n_jets_fid" , 							&tn_jets_fid ,      					"n_jets_fid/I"   				);		
	tree->Branch("leading_jet_pt_fid" , 					&tleading_jet_pt_fid ,      			"leading_jet_pt_fid/F"   		);	
	tree->Branch("n_jets_truth_fid" , 						&tn_jets_truth_fid ,      				"n_jets_truth_fid/I"   			);		
	tree->Branch("leading_jet_pt_truth_fid" , 				&tleading_jet_pt_truth_fid ,      		"leading_jet_pt_truth_fid/F"   	);	


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
	
	tree->Branch("valTrackIso" , 							&tValTrackIso ,           				"tValTrackIso[4]/F"  	);
	tree->Branch("valCaloIso" , 							&tValCaloIso ,           				"tValCaloIso[4]/F"  	);
	tree->Branch("valD0Sign" , 								&tValD0Sign ,           				"tValD0Sign[4]/F"  		);
	tree->Branch("trackIso" , 								&tTrackIso ,           					"tTrackIso[4]/O" 	 	);
	tree->Branch("caloIso" , 								&tCaloIso ,           					"tCaloIso[4]/O"  		);
	tree->Branch("d0Sign" , 								&tD0Sign ,           					"tD0Sign[4]/O"  		);
	tree->Branch("ElID" , 									&tElID ,        	   					"tElID[4]/O"  		);

	tree->Branch("flagQuad" , 								&tflagQuad ,           					"flagQuad/I"  			);	
<<<<<<< HEAD
=======

	tree->Branch("BCHCutMedium" , 							&tBCHCutMedium ,      					"BCHCutMedium/I"   );		
	tree->Branch("BCHCutTight" , 							&tBCHCutTight ,      					"BCHCutTight/I"   );		
	
>>>>>>> bc7b9ddaf72f0a41dfe1bb5d9068cc4b03444c0d
}

void OutputTree::fillTree(D3PDReader::Event *event, QuadLepton * higgs, Int_t type, Bool_t isMC)
{
	// Filling the variables
	if(isMC) tRun = event->eventinfo.mc_channel_number();
	else tRun = event->eventinfo.RunNumber();
	tEvent = event->eventinfo.EventNumber();;
	tlbn = event->eventinfo.lbn();

	tm4l_unconstrained 		= higgs->getMass()/1000;
	tm4lerr_unconstrained 	= higgs->getMassErr()/1000;
	tmZ1_unconstrained 		= higgs->getZ1()->get4Momentum()->M()/1000;
	tmZ2_unconstrained 		= higgs->getZ2()->get4Momentum()->M()/1000;

	tm4l_fsr 	= higgs->getMassFSR()/1000;
	tmZ1_fsr 	= higgs->getZ1MassFSR()/1000;
	tmZ2_fsr 	= higgs->getZ2MassFSR()/1000;
	tm4lerr_fsr = higgs->getMassErrFSR()/1000;
	tfsrType	= higgs->fsrType;

	tm4l_constrained 		= higgs->getMassZMassCons()/1000;
	tm4lerr_constrained 	= higgs->getMassErrZMassCons()/1000;
	tmZ1_constrained 		= higgs->getZ1MassZMassCons()/1000;
	tmZ2_constrained 		= higgs->getZ2MassZMassCons()/1000;

	tm4l_unconstrained_ID 		= higgs->massID/1000;
	tm4lerr_unconstrained_ID 	= higgs->massErrID/1000;
	tm4l_constrained_ID 		= higgs->massZMassConsID/1000;
	tm4lerr_constrained_ID	 	= higgs->massErrZmassConsID/1000;
	tm4l_unconstrained_MS 		= higgs->massMS/1000;
	tm4lerr_unconstrained_MS 	= higgs->massErrMS/1000;
	tm4l_constrained_MS 		= higgs->massZMassConsMS/1000;
	tm4lerr_constrained_MS	 	= higgs->massErrZmassConsMS/1000;

	tweight 		= higgs->weight;
	tweight_corr 	= higgs->weight_corr;
	tweight_lumi 	= higgs->weight_lumi;

	tpt4l_unconstrained 	= higgs->sum_unconstrained.Pt()/1000;
	ty4l_unconstrained 		= higgs->sum_unconstrained.Rapidity();
	teta4l_unconstrained 	= higgs->sum_unconstrained.Eta();
	tpt4l_fsr 				= higgs->sum_fsr.Pt()/1000;
	ty4l_fsr 				= higgs->sum_fsr.Rapidity();
	teta4l_fsr 				= higgs->sum_fsr.Eta();  
	tpt4l_constrained 		= higgs->sum_constrained.Pt()/1000;
	ty4l_constrained 		= higgs->sum_constrained.Rapidity();
	teta4l_constrained 		= higgs->sum_constrained.Eta();

	if(isMC && higgs->truthVec.Pt() != 0)
	{
		tpt4l_truth_born 		= higgs->truthVec.Pt()/1000;
		ty4l_truth_born 		= higgs->truthVec.Rapidity();
		teta4l_truth_born 		= higgs->truthVec.Eta();
<<<<<<< HEAD
=======
		tphi4l_truth_born 		= higgs->truthVec.Phi();
>>>>>>> bc7b9ddaf72f0a41dfe1bb5d9068cc4b03444c0d
		tm4l_truth_born 		= higgs->truthVec.M()/1000;
	}
	tZ1_lepplus_pt		= higgs->getZ1()->getLepPlus()->get4Momentum()->Pt()/1000;
	tZ1_lepminus_pt 	= higgs->getZ1()->getLepNeg()->get4Momentum()->Pt()/1000;
	tZ2_lepplus_pt 		= higgs->getZ2()->getLepPlus()->get4Momentum()->Pt()/1000;
	tZ2_lepminus_pt 	= higgs->getZ2()->getLepNeg()->get4Momentum()->Pt()/1000;

	tZ1_lepplus_eta 	= higgs->getZ1()->getLepPlus()->get4Momentum()->Eta();	
	tZ1_lepminus_eta 	= higgs->getZ1()->getLepNeg()->get4Momentum()->Eta();
	tZ2_lepplus_eta 	= higgs->getZ2()->getLepPlus()->get4Momentum()->Eta();
	tZ2_lepminus_eta 	= higgs->getZ2()->getLepNeg()->get4Momentum()->Eta();

	tZ1_lepplus_phi 	= higgs->getZ1()->getLepPlus()->get4Momentum()->Phi();
	tZ1_lepminus_phi 	= higgs->getZ1()->getLepNeg()->get4Momentum()->Phi();
	tZ2_lepplus_phi 	= higgs->getZ2()->getLepPlus()->get4Momentum()->Phi();
	tZ2_lepminus_phi 	= higgs->getZ2()->getLepNeg()->get4Momentum()->Phi();

	tZ1_lepplus_m 		= higgs->getZ1()->getLepPlus()->lepMass/1000;
	tZ1_lepminus_m 		= higgs->getZ1()->getLepNeg()->lepMass/1000;
	tZ2_lepplus_m 		= higgs->getZ2()->getLepPlus()->lepMass/1000;
	tZ2_lepminus_m 		= higgs->getZ2()->getLepNeg()->lepMass/1000;

	tZ1_lepplus_id 		= (Int_t) higgs->getZ1()->getLepPlus()->lepID;
	tZ1_lepminus_id 	= (Int_t) higgs->getZ1()->getLepNeg()->lepID;
	tZ2_lepplus_id 		= (Int_t) higgs->getZ2()->getLepPlus()->lepID;
	tZ2_lepminus_id 	= (Int_t) higgs->getZ2()->getLepNeg()->lepID;
<<<<<<< HEAD
	
=======

	tZ1_lepplus_turthParent 	= (Int_t) higgs->getZ1()->getLepPlus()->truthParentType;
	tZ1_lepminus_turthParent 	= (Int_t) higgs->getZ1()->getLepNeg()->truthParentType;
	tZ2_lepplus_turthParent 	= (Int_t) higgs->getZ2()->getLepPlus()->truthParentType;
	tZ2_lepminus_turthParent 	= (Int_t) higgs->getZ2()->getLepNeg()->truthParentType;

>>>>>>> bc7b9ddaf72f0a41dfe1bb5d9068cc4b03444c0d
	tZ1_lepplus_cov_mom  	= higgs->getZ1()->getLepPlus()->covMomErr;
	tZ1_lepminus_cov_mom  	= higgs->getZ1()->getLepNeg()->covMomErr;
	tZ2_lepplus_cov_mom  	= higgs->getZ2()->getLepPlus()->covMomErr;
	tZ2_lepminus_cov_mom  	= higgs->getZ2()->getLepNeg()->covMomErr;
	

	tZ1_lepplus_pt_uncorr  		= higgs->getZ1()->getLepPlus()->cb_pt_unsmeared/1000;
	tZ1_lepminus_pt_uncorr  	= higgs->getZ1()->getLepNeg()->cb_pt_unsmeared/1000;
	tZ2_lepplus_pt_uncorr 	 	= higgs->getZ2()->getLepPlus()->cb_pt_unsmeared/1000;
	tZ2_lepminus_pt_uncorr  	= higgs->getZ2()->getLepNeg()->cb_pt_unsmeared/1000;
	
	tZ1_lepplus_pt_uncorr_ID  	= higgs->getZ1()->getLepPlus()->id_pt_unsmeared/1000;
	tZ1_lepminus_pt_uncorr_ID  	= higgs->getZ1()->getLepNeg()->id_pt_unsmeared/1000;
	tZ2_lepplus_pt_uncorr_ID 	= higgs->getZ2()->getLepPlus()->id_pt_unsmeared/1000;
	tZ2_lepminus_pt_uncorr_ID  	= higgs->getZ2()->getLepNeg()->id_pt_unsmeared/1000;
	
	tZ1_lepplus_pt_uncorr_MS 	= higgs->getZ1()->getLepPlus()->me_pt_unsmeared/1000;
	tZ1_lepminus_pt_uncorr_MS  	= higgs->getZ1()->getLepNeg()->me_pt_unsmeared/1000;
	tZ2_lepplus_pt_uncorr_MS  	= higgs->getZ2()->getLepPlus()->me_pt_unsmeared/1000;
	tZ2_lepminus_pt_uncorr_MS  	= higgs->getZ2()->getLepNeg()->me_pt_unsmeared/1000;


	if(higgs->getZ1()->getLepPlus()->m_momentumTruthRecoBare.Pt() && isMC)
	{
		tZ1_lepplus_pt_truth_bare 	= higgs->getZ1()->getLepPlus()->m_momentumTruthRecoBare.Pt()/1000;
		tZ1_lepplus_eta_truth_bare 	= higgs->getZ1()->getLepPlus()->m_momentumTruthRecoBare.Eta();
		tZ1_lepplus_phi_truth_bare 	= higgs->getZ1()->getLepPlus()->m_momentumTruthRecoBare.Phi();
		tZ1_lepplus_m_truth_bare 	= higgs->getZ1()->getLepPlus()->m_momentumTruthRecoBare.M()/1000;
	}
	if(higgs->getZ1()->getLepNeg()->m_momentumTruthRecoBare.Pt() && isMC)
	{
		tZ1_lepminus_pt_truth_bare 	= higgs->getZ1()->getLepNeg()->m_momentumTruthRecoBare.Pt()/1000;
		tZ1_lepminus_eta_truth_bare = higgs->getZ1()->getLepNeg()->m_momentumTruthRecoBare.Eta();
		tZ1_lepminus_phi_truth_bare = higgs->getZ1()->getLepNeg()->m_momentumTruthRecoBare.Phi();
		tZ1_lepminus_m_truth_bare 	= higgs->getZ1()->getLepNeg()->m_momentumTruthRecoBare.M()/1000;
	}
	if(higgs->getZ2()->getLepPlus()->m_momentumTruthRecoBare.Pt() && isMC)
	{
		tZ2_lepplus_pt_truth_bare 	= higgs->getZ2()->getLepPlus()->m_momentumTruthRecoBare.Pt()/1000;
		tZ2_lepplus_eta_truth_bare 	= higgs->getZ2()->getLepPlus()->m_momentumTruthRecoBare.Eta();
		tZ2_lepplus_phi_truth_bare 	= higgs->getZ2()->getLepPlus()->m_momentumTruthRecoBare.Phi();
		tZ2_lepplus_m_truth_bare 	= higgs->getZ2()->getLepPlus()->m_momentumTruthRecoBare.M()/1000;
	}
	if(higgs->getZ2()->getLepNeg()->m_momentumTruthRecoBare.Pt() && isMC)
	{
		tZ2_lepminus_pt_truth_bare 	= higgs->getZ2()->getLepNeg()->m_momentumTruthRecoBare.Pt()/1000;
		tZ2_lepminus_eta_truth_bare = higgs->getZ2()->getLepNeg()->m_momentumTruthRecoBare.Eta();
		tZ2_lepminus_phi_truth_bare = higgs->getZ2()->getLepNeg()->m_momentumTruthRecoBare.Phi();
		tZ2_lepminus_m_truth_bare 	= higgs->getZ2()->getLepNeg()->m_momentumTruthRecoBare.M()/1000;
	}
  
	if(higgs->getZ1()->getLepPlus()->m_momentumTruthTrueBorn.Pt() && isMC)
	{
		tZ1_lepplus_pt_truth 	= higgs->getZ1()->getLepPlus()->m_momentumTruthTrueBorn.Pt()/1000;
		tZ1_lepplus_eta_truth 	= higgs->getZ1()->getLepPlus()->m_momentumTruthTrueBorn.Eta();
		tZ1_lepplus_phi_truth 	= higgs->getZ1()->getLepPlus()->m_momentumTruthTrueBorn.Phi();
		tZ1_lepplus_m_truth 	= higgs->getZ1()->getLepPlus()->m_momentumTruthTrueBorn.M()/1000;
	}
	if(higgs->getZ1()->getLepNeg()->m_momentumTruthTrueBorn.Pt() && isMC)
	{
		tZ1_lepminus_pt_truth	= higgs->getZ1()->getLepNeg()->m_momentumTruthTrueBorn.Pt()/1000;
		tZ1_lepminus_eta_truth 	= higgs->getZ1()->getLepNeg()->m_momentumTruthTrueBorn.Eta();
		tZ1_lepminus_phi_truth 	= higgs->getZ1()->getLepNeg()->m_momentumTruthTrueBorn.Phi();
		tZ1_lepminus_m_truth 	= higgs->getZ1()->getLepNeg()->m_momentumTruthTrueBorn.M()/1000;
	}
	if(higgs->getZ2()->getLepPlus()->m_momentumTruthTrueBorn.Pt() && isMC)
	{
		tZ2_lepplus_pt_truth 	= higgs->getZ2()->getLepPlus()->m_momentumTruthTrueBorn.Pt()/1000;
		tZ2_lepplus_eta_truth 	= higgs->getZ2()->getLepPlus()->m_momentumTruthTrueBorn.Eta();
		tZ2_lepplus_phi_truth 	= higgs->getZ2()->getLepPlus()->m_momentumTruthTrueBorn.Phi();
		tZ2_lepplus_m_truth 	= higgs->getZ2()->getLepPlus()->m_momentumTruthTrueBorn.M()/1000;
	}
	if(higgs->getZ2()->getLepNeg()->m_momentumTruthTrueBorn.Pt() && isMC)
	{
		tZ2_lepminus_pt_truth 	= higgs->getZ2()->getLepNeg()->m_momentumTruthTrueBorn.Pt()/1000;
		tZ2_lepminus_eta_truth 	= higgs->getZ2()->getLepNeg()->m_momentumTruthTrueBorn.Eta();
		tZ2_lepminus_phi_truth 	= higgs->getZ2()->getLepNeg()->m_momentumTruthTrueBorn.Phi();
		tZ2_lepminus_m_truth	= higgs->getZ2()->getLepNeg()->m_momentumTruthTrueBorn.M()/1000;
	}

	tn_jets		 				= higgs->n_jets;
	tdijet_invmass		 		= higgs->dijet_invmass;
	tdijet_deltaeta 			= higgs->dijet_deltaeta;
	if(higgs->leadingJet)
	{
		tleading_jet_pt			= higgs->leadingJet->get4Momentum()->Pt()/1000;
		tleading_jet_eta		= higgs->leadingJet->get4Momentum()->Eta();
		tleading_jet_phi		= higgs->leadingJet->get4Momentum()->Phi();
		tleading_jet_m			= higgs->leadingJet->get4Momentum()->M()/1000;
		tleading_jet_width		= higgs->leadingJet->GetJets()->WIDTH();
		tleading_jet_nTrk		= higgs->leadingJet->GetJets()->nTrk();
<<<<<<< HEAD
=======
		//cout<<"Output: Leading_jet_m "<<tleading_jet_m<<endl;
>>>>>>> bc7b9ddaf72f0a41dfe1bb5d9068cc4b03444c0d
	}
	if(higgs->subLeadingJet)
	{
		tsubleading_jet_pt		= higgs->subLeadingJet->get4Momentum()->Pt()/1000;
		tsubleading_jet_eta		= higgs->subLeadingJet->get4Momentum()->Eta();
		tsubleading_jet_phi		= higgs->subLeadingJet->get4Momentum()->Phi();
		tsubleading_jet_m		= higgs->subLeadingJet->get4Momentum()->M()/1000;
		tsubleading_jet_width	= higgs->subLeadingJet->GetJets()->WIDTH();
		tsubleading_jet_nTrk	= higgs->subLeadingJet->GetJets()->nTrk();
	}
	if(higgs->thirdJet)
	{
		tthird_jet_pt			= higgs->thirdJet->get4Momentum()->Pt()/1000;
		tthird_jet_eta			= higgs->thirdJet->get4Momentum()->Eta();
		tthird_jet_phi			= higgs->thirdJet->get4Momentum()->Phi();
		tthird_jet_m			= higgs->thirdJet->get4Momentum()->M()/1000;
		tthird_jet_width		= higgs->thirdJet->GetJets()->WIDTH();
		tthird_jet_nTrk			= higgs->thirdJet->GetJets()->nTrk();
	}
	//tleading_jet_pt 			= higgs->leading_jet_pt;
	//tleading_jet_eta 			= higgs->leading_jet_eta;
	//tsubleading_jet_pt 		= higgs->subleading_jet_pt;
	tBDT_discriminant_VBF		= higgs->BDT_discriminant_VBF;
	tBDT_discriminant_HadVH 	= higgs->BDT_discriminant_HadVH;
	tn_jets_truth_bare 			= higgs->n_jets_truth_bare;
	tleading_jet_pt_truth_bare 	= higgs->leading_jet_pt_truth_bare;
	
	tn_jets_fid		 			= higgs->n_jets_fid;
	tleading_jet_pt_fid		 	= higgs->leading_jet_pt_fid;
	tn_jets_truth_fid		 	= higgs->n_jets_truth_fid;
	tleading_jet_pt_truth_fid	= higgs->leading_jet_pt_truth_fid;

	tcthstr 	= higgs->cthstr;
	tphi1 		= higgs->phi1;
	tcth1		= higgs->cth1;
	tcth2 		= higgs->cth2;
	tphi 		= higgs->phi;
	thiggspt 	= higgs->get4Momentum()->Pt()/1000;
	thiggseta 	= higgs->get4Momentum()->Eta();

	tzVertWeight 	= higgs->zVertWeight;
	tpileupWeight 	= higgs->pileupWeight;
	tggFWeight 		= higgs->ggFWeight;
	tJHUWeight 		= higgs->JHUWeight;	
	tttBarVeto 		= higgs->ttBarVeto;
	tzzOverlapVeto 	= higgs->zzOverlapVeto;
	tzBBVeto 		= higgs->zBBVeto;
	ttrigEff 		= higgs->trigEff;
	tlepEff 		= higgs->lepEff;
	tmcEventWeight 	= higgs->mcEventWeight;
	tcrossSection 	= higgs->crossSection;
	tbranchRatio 	= higgs->branchRatio;
	tlumi			= higgs->lumi;
	tflagQuad 		= higgs->flagQuad;

	// Truth
	tm4l_truth_matched_born 	= higgs->m4l_truth_born;
	tmZ1_truth_matched_born 	= higgs->mZ1_truth_born;
	tmZ2_truth_matched_born		= higgs->mZ2_truth_born;
	tm4l_truth_matched_bare		= higgs->m4l_truth_bare;
	tmZ1_truth_matched_bare 	= higgs->mZ1_truth_bare;
	tmZ2_truth_matched_bare 	= higgs->mZ2_truth_bare;
	tm4l_truth_matched_dressed 	= higgs->m4l_truth_dressed;
	tmZ1_truth_matched_dressed 	= higgs->mZ1_truth_dressed;
	tmZ2_truth_matched_dressed 	= higgs->mZ2_truth_dressed;

	// BDT 
	tKD_discriminant		= higgs->KD_discriminant;
	tBDT_discriminant		= higgs->BDT_discriminant;
	tBDTGuass_discriminant	= higgs->BDTGuass_discriminant;
	tptSysupFac				= higgs->ptSysupFac;
	tptSysdownFac			= higgs->ptSysdownFac;

	// For Category
<<<<<<< HEAD
	tmissing_et = event->MET_RefFinal.et();
=======
	tmissing_et = event->MET_RefFinal.et()/1000;
>>>>>>> bc7b9ddaf72f0a41dfe1bb5d9068cc4b03444c0d
	if(higgs->leadingExtraLep)
	{
		tleading_additional_lepton_pt 		= higgs->leadingExtraLep->get4Momentum()->Pt()/1000;
		tleading_additional_lepton_eta 		= higgs->leadingExtraLep->get4Momentum()->Eta();
		tleading_additional_lepton_phi 		= higgs->leadingExtraLep->get4Momentum()->Phi();
		tleading_additional_lepton_m 		= higgs->leadingExtraLep->get4Momentum()->M()/1000;
		tleading_additional_lepton_type		= higgs->leadingExtraLep->lepType;
		tleading_additional_lepton_type_truth_matched_bare = higgs->leadingExtraLep->truthParentType;
	}
	if(higgs->subleadingExtraLep)
	{
		tsubleading_additional_lepton_pt 	= higgs->subleadingExtraLep->get4Momentum()->Pt()/1000;
		tsubleading_additional_lepton_eta 	= higgs->subleadingExtraLep->get4Momentum()->Eta();
		tsubleading_additional_lepton_phi 	= higgs->subleadingExtraLep->get4Momentum()->Phi();
		tsubleading_additional_lepton_m 	= higgs->subleadingExtraLep->get4Momentum()->M()/1000;
		tsubleading_additional_lepton_type	= higgs->subleadingExtraLep->lepType;
		tsubleading_additional_lepton_type_truth_matched_bare = higgs->subleadingExtraLep->truthParentType;
	}

	// NPV
	tnpv	= higgs->npv;

	// Calibration
	tcalib	= higgs->calib;

	// eventType
	teventType = type;

	// for CR
	for(Int_t i = 0; i < 4; i++)
	{
		tValTrackIso[i] = higgs->valTrackIso[i];
		tValCaloIso[i] 	= higgs->valCaloIso[i];
		tValD0Sign[i] 	= higgs->valD0Sig[i];

		tTrackIso[i]	= higgs->trackIso[i];
		tCaloIso[i] 	= higgs->caloIso[i];
		tD0Sign[i] 		= higgs->d0Sig[i];
		tElID[i]		= higgs->looseElectron[i];
	}

<<<<<<< HEAD
=======
	tBCHCutMedium = higgs->BCHCutMedium;
	tBCHCutTight = higgs->BCHCutTight;

	//cout<<"Ouput BCH Medium "<<tBCHCutMedium<<endl;
	//cout<<"Ouput BCH tight "<<tBCHCutTight<<endl;

>>>>>>> bc7b9ddaf72f0a41dfe1bb5d9068cc4b03444c0d
	// Filling the tree
	if(type == analysisType::Mu4) mu4Tree->Fill();
	else if(type == analysisType::El4) el4Tree->Fill();
	else if(type == analysisType::Mu2El2) mu2el2Tree->Fill();
	else if(type == analysisType::El2Mu2) el2mu2Tree->Fill();
	else {cout<<"OutputTree: fillTree: Analysis type not recognized"<<endl;}
	// Production Category
	if(type == analysisType::Mu4 && higgs->prodChannel == productionChannel::ggF) mu4ggFTree->Fill();
	else if(type == analysisType::El4 && higgs->prodChannel == productionChannel::ggF) el4ggFTree->Fill();
	else if(type == analysisType::Mu2El2 && higgs->prodChannel == productionChannel::ggF) mu2el2ggFTree->Fill();
	else if(type == analysisType::El2Mu2 && higgs->prodChannel == productionChannel::ggF) el2mu2ggFTree->Fill();
	else if(higgs->prodChannel == productionChannel::VBF) VBFTree->Fill();
	else if(higgs->prodChannel == productionChannel::VHLep) VHLepTree->Fill();
	else if(higgs->prodChannel == productionChannel::VHHad) VHHadTree->Fill();	
	else if(higgs->prodChannel == productionChannel::VH) VHTree->Fill();	
	else {cout<<"OutputTree: fillTree: Analysis type not recognized Prod Category"<<endl;}

	// Clearing the vars for the next one
	clearVars();
}

void OutputTree::saveTrees(TString filePath, TH1F* countingHist, TString sampleName)
{
	TFile *output = new TFile (filePath, "RECREATE");

	output->cd();
	// Counting Hist
	countingHist->Write();
	// Saving the trees
	mu4Tree->Write();
	mu2el2Tree->Write();
	el2mu2Tree->Write();
	el4Tree->Write();
	
	mu4ggFTree->Write();
	mu2el2ggFTree->Write();
	el2mu2ggFTree->Write();
	el4ggFTree->Write();
	
	VBFTree->Write();
	VHHadTree->Write();
	VHLepTree->Write();
	VHTree->Write();
	// Output Folder with the sample name for later checking
	output->mkdir(sampleName);
	cout<<"Sample Name: "<<sampleName<<endl;

	output->Close();


}


