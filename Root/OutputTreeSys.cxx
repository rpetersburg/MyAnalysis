#include <stdlib.h>
#include <string>
#include "MyAnalysis/OutputTreeSys.h"


using namespace std;
////////////////////////////////////////////////////////////////////////////////////////
//				Constructor and Festructor
////////////////////////////////////////////////////////////////////////////////////////
OutputTreeSys::OutputTreeSys()
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

OutputTreeSys::~OutputTreeSys()
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

void OutputTreeSys::clearVars()
{
	tRun = -999;
	tEvent = -999;
	tlbn = -999;

	tm4l_unconstrained = -999;
	tm4lerr_unconstrained = -999;
	tmZ1_unconstrained = -999;
	tmZ2_unconstrained = -999;

	tm4l_constrained = -999;
	tm4lerr_constrained = -999;
	tmZ1_constrained = -999;
	tmZ2_constrained = -999;

	tweight = -999;
	tweight_corr = -999;
	tweight_lumi = -999;


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
	
	tZ1_lepplus_pt_uncorr_CB = -999;
	tZ1_lepminus_pt_uncorr_CB = -999;
	tZ2_lepplus_pt_uncorr_CB = -999;
	tZ2_lepminus_pt_uncorr_CB = -999;

	tZ1_lepplus_pt_uncorr_ID = -999;
	tZ1_lepminus_pt_uncorr_ID = -999;
	tZ2_lepplus_pt_uncorr_ID = -999;
	tZ2_lepminus_pt_uncorr_ID = -999;

	tZ1_lepplus_pt_uncorr_MS = -999;
	tZ1_lepminus_pt_uncorr_MS = -999;
	tZ2_lepplus_pt_uncorr_MS = -999;
	tZ2_lepminus_pt_uncorr_MS = -999;

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

	tflagQuad = -999;	
}

void OutputTreeSys::bookTree(TTree *tree)
{
	//cout<<"Booking Tree: "<<tree->getName() <<endl;

	tree->Branch("run" , 									&tRun , 								"run/I"					);		
	tree->Branch("event" , 									&tEvent ,                   			"event/I" 				);	
	tree->Branch("lbn" , 									&tlbn ,                     			"lbn/I" 				);	
	
	tree->Branch("m4l_unconstrained" , 						&tm4l_unconstrained ,       			"m4l_unconstrained/F"   );
	tree->Branch("m4lerr_unconstrained" , 					&tm4lerr_unconstrained ,    			"m4lerr_unconstrained/F"); 	
	tree->Branch("mZ1_unconstrained" , 						&tmZ1_unconstrained ,       			"mZ1_unconstrained/F"	);	
	tree->Branch("mZ2_unconstrained" , 						&tmZ2_unconstrained ,       			"mZ2_unconstrained/F" 	);	
	
	tree->Branch("m4l_constrained" , 						&tm4l_constrained ,         			"m4l_constrained/F"     );
	tree->Branch("m4lerr_constrained" ,						&tm4lerr_constrained ,      			"m4lerr_constrained/F" 	);	
	tree->Branch("mZ1_constrained" , 						&tmZ1_constrained ,         			"mZ1_constrained/F"  	);	
	tree->Branch("mZ2_constrained" , 						&tmZ2_constrained ,         			"mZ2_constrained/F"  	);	

	tree->Branch("weight" , 								&tweight ,                  			"weight/F"  			);	
	tree->Branch("weight_corr" , 							&tweight_corr ,             			"weight_corr/F"  		);	
	tree->Branch("weight_lumi" , 							&tweight_lumi ,             			"weight_lumi/F"  		);	

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
	tree->Branch("Z1_lepplus_E" , 							&tZ1_lepplus_E ,            			"Z1_lepplus_E/F"  		);	
	tree->Branch("Z1_lepminus_E" , 							&tZ1_lepminus_E ,           			"Z1_lepminus_E/F"  		);	
	tree->Branch("Z2_lepplus_E" , 							&tZ2_lepplus_E ,            			"Z2_lepplus_E/F"  		);	
	tree->Branch("Z2_lepminus_E" , 							&tZ2_lepminus_E ,           			"Z2_lepminus_E/F"  		);
	
	tree->Branch("Z1_lepplus_cov_mom" , 					&tZ1_lepplus_cov_mom ,      			"Z1_lepplus_cov_mom/F"   );	
	tree->Branch("Z1_lepminus_cov_mom" , 					&tZ1_lepminus_cov_mom ,     			"Z1_lepminus_cov_mom/F"  );	
	tree->Branch("Z2_lepplus_cov_mom" , 					&tZ2_lepplus_cov_mom ,      			"Z2_lepplus_cov_mom/F"   );	
	tree->Branch("Z2_lepminus_cov_mom" , 					&tZ2_lepminus_cov_mom ,     			"Z2_lepminus_cov_mom/F"  );	

	tree->Branch("Z1_lepplus_id" , 							&tZ1_lepplus_id ,            			"Z1_lepplus_id/I"  		);	
	tree->Branch("Z1_lepminus_id" , 						&tZ1_lepminus_id ,           			"Z1_lepminus_id/I"  	);	
	tree->Branch("Z2_lepplus_id" , 							&tZ2_lepplus_id ,            			"Z2_lepplus_id/I"  		);	
	tree->Branch("Z2_lepminus_id" , 						&tZ2_lepminus_id ,           			"Z2_lepminus_id/I"  	);

	tree->Branch("Z1_lepplus_pt_uncorr_CB" , 				&tZ1_lepplus_pt_uncorr_CB ,           	"Z1_lepplus_pt_uncorr_CB/F"  	);	
	tree->Branch("Z1_lepminus_pt_uncorr_CB" , 				&tZ1_lepminus_pt_uncorr_CB ,          	"Z1_lepminus_pt_uncorr_CB/F"  	);	
	tree->Branch("Z2_lepplus_pt_uncorr_CB" , 				&tZ2_lepplus_pt_uncorr_CB ,           	"Z2_lepplus_pt_uncorr_CB/F" 	);	
	tree->Branch("Z2_lepminus_pt_uncorr_CB" , 				&tZ2_lepminus_pt_uncorr_CB ,          	"Z2_lepminus_pt_uncorr_CB/F"  	);

	tree->Branch("Z1_lepplus_pt_uncorr_ID" , 				&tZ1_lepplus_pt_uncorr_ID ,           	"Z1_lepplus_pt_uncorr_ID/F"  	);	
	tree->Branch("Z1_lepminus_pt_uncorr_ID" , 				&tZ1_lepminus_pt_uncorr_ID ,          	"Z1_lepminus_pt_uncorr_ID/F"  	);	
	tree->Branch("Z2_lepplus_pt_uncorr_ID" , 				&tZ2_lepplus_pt_uncorr_ID ,           	"Z2_lepplus_pt_uncorr_ID/F" 	);	
	tree->Branch("Z2_lepminus_pt_uncorr_ID" , 				&tZ2_lepminus_pt_uncorr_ID ,          	"Z2_lepminus_pt_uncorr_ID/F"  	);
	
	tree->Branch("Z1_lepplus_pt_uncorr_MS" , 				&tZ1_lepplus_pt_uncorr_MS ,           	"Z1_lepplus_pt_uncorr_MS/F"  	);	
	tree->Branch("Z1_lepminus_pt_uncorr_MS" , 				&tZ1_lepminus_pt_uncorr_MS ,          	"Z1_lepminus_pt_uncorr_MS/F"  	);	
	tree->Branch("Z2_lepplus_pt_uncorr_MS" , 				&tZ2_lepplus_pt_uncorr_MS ,           	"Z2_lepplus_pt_uncorr_MS/F" 	);	
	tree->Branch("Z2_lepminus_pt_uncorr_MS" , 				&tZ2_lepminus_pt_uncorr_MS ,          	"Z2_lepminus_pt_uncorr_MS/F"  	);

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

	tree->Branch("Lep_cl_eta" ,								&Lep_cl_eta,							"Lep_cl_eta[4]/F"					);
	tree->Branch("Lep_E_ZeeStatUp" ,						&Lep_E_ZeeStatUp,						"Lep_E_ZeeStatUp[4]/F"				);
	tree->Branch("Lep_E_ZeeStatDown" ,						&Lep_E_ZeeStatDown,						"Lep_E_ZeeStatDown[4]/F"			);
	tree->Branch("Lep_E_ZeeSystUp" ,						&Lep_E_ZeeSystUp,						"Lep_E_ZeeSystUp[4]/F"				);
	tree->Branch("Lep_E_ZeeSystDown" ,						&Lep_E_ZeeSystDown,						"Lep_E_ZeeSystDown[4]/F"			);	
	tree->Branch("Lep_E_ZeeAllUp" ,							&Lep_E_ZeeAllUp,						"Lep_E_ZeeAllUp[4]/F"				);
	tree->Branch("Lep_E_ZeeAllDown" ,						&Lep_E_ZeeAllDown,						"Lep_E_ZeeAllDown[4]/F"				);
	tree->Branch("Lep_E_PSUp" ,								&Lep_E_PSUp,							"Lep_E_PSUp[4]/F"					);
	tree->Branch("Lep_E_PSDown" ,							&Lep_E_PSDown,							"Lep_E_PSDown[4]/F"					);
	tree->Branch("Lep_E_S12Up" ,							&Lep_E_S12Up,							"Lep_E_S12Up[4]/F"					);
	tree->Branch("Lep_E_S12Down" ,							&Lep_E_S12Down,							"Lep_E_S12Down[4]/F"				);
	tree->Branch("Lep_E_MatIDUp" ,							&Lep_E_MatIDUp,							"Lep_E_MatIDUp[4]/F"				);
	tree->Branch("Lep_E_MatIDDown" ,						&Lep_E_MatIDDown,						"Lep_E_MatIDDown[4]/F"				);
	tree->Branch("Lep_E_MatCryoUp" ,						&Lep_E_MatCryoUp,						"Lep_E_MatCryoUp[4]/F"				);
	tree->Branch("Lep_E_MatCryoDown" ,						&Lep_E_MatCryoDown,						"Lep_E_MatCryoDown[4]/F"			);
	tree->Branch("Lep_E_MatCaloUp" ,						&Lep_E_MatCaloUp,						"Lep_E_MatCaloUp[4]/F"				);
	tree->Branch("Lep_E_MatCaloDown" ,						&Lep_E_MatCaloDown,						"Lep_E_MatCaloDown[4]/F"			);
	tree->Branch("Lep_E_LArCalibUp" ,						&Lep_E_LArCalibUp,						"Lep_E_LArCalibUp[4]/F"				);
	tree->Branch("Lep_E_LArCalibDown" ,						&Lep_E_LArCalibDown,					"Lep_E_LArCalibDown[4]/F"			);
	tree->Branch("Lep_E_LArUnconvCalibUp" ,					&Lep_E_LArUnconvCalibUp,				"Lep_E_LArUnconvCalibUp[4]/F"		);
	tree->Branch("Lep_E_LArUnconvCalibDown" ,				&Lep_E_LArUnconvCalibDown,				"Lep_E_LArUnconvCalibDown[4]/F"		);
	tree->Branch("Lep_E_LArElecUnconvUp" ,					&Lep_E_LArElecUnconvUp,					"Lep_E_LArElecUnconvUp[4]/F"		);
	tree->Branch("Lep_E_LArElecUnconvDown" ,				&Lep_E_LArElecUnconvDown,				"Lep_E_LArElecUnconvDown[4]/F"		);
	tree->Branch("Lep_E_LArElecCalibUp" ,					&Lep_E_LArElecCalibUp,					"Lep_E_LArElecCalibUp[4]/F"			);
	tree->Branch("Lep_E_LArElecCalibDown" ,					&Lep_E_LArElecCalibDown,				"Lep_E_LArElecCalibDown[4]/F"		);
	tree->Branch("Lep_E_GainUp" ,							&Lep_E_GainUp,							"Lep_E_GainUp[4]/F"					);
	tree->Branch("Lep_E_GainDown" ,							&Lep_E_GainDown,						"Lep_E_GainDown[4]/F"				);
	tree->Branch("Lep_E_G4Up" ,								&Lep_E_G4Up,							"Lep_E_G4Up[4]/F"					);
	tree->Branch("Lep_E_G4Down" ,							&Lep_E_G4Down,							"Lep_E_G4Down[4]/F"					);
	tree->Branch("Lep_E_MomentumUp" ,						&Lep_E_MomentumUp,						"Lep_E_MomentumUp[4]/F"				);
	tree->Branch("Lep_E_MomentumDown" ,						&Lep_E_MomentumDown,					"Lep_E_MomentumDown[4]/F"			);
	tree->Branch("Lep_E_ZSmearingUp" ,						&Lep_E_ZSmearingUp,						"Lep_E_ZSmearingUp[4]/F"			);
	tree->Branch("Lep_E_ZSmearingDown" ,					&Lep_E_ZSmearingDown,					"Lep_E_ZSmearingDown[4]/F"			);
	tree->Branch("Lep_E_SamplingTermUp" ,					&Lep_E_SamplingTermUp,					"Lep_E_SamplingTermUp[4]/F"			);
	tree->Branch("Lep_E_SamplingTermDown" ,					&Lep_E_SamplingTermDown,				"Lep_E_SamplingTermDown[4]/F"		);
	tree->Branch("Lep_E_MaterialIDUp" ,						&Lep_E_MaterialIDUp,					"Lep_E_MaterialIDUp[4]/F"			);
	tree->Branch("Lep_E_MaterialIDDown" ,					&Lep_E_MaterialIDDown,					"Lep_E_MaterialIDDown[4]/F"			);
	tree->Branch("Lep_E_MaterialCaloUp" ,					&Lep_E_MaterialCaloUp,					"Lep_E_MaterialCaloUp[4]/F"			);
	tree->Branch("Lep_E_MaterialCaloDown" ,					&Lep_E_MaterialCaloDown,				"Lep_E_MaterialCaloDown[4]/F"		);
	tree->Branch("Lep_E_MaterialGapUp" ,					&Lep_E_MaterialGapUp,					"Lep_E_MaterialGapUp[4]/F"			);
	tree->Branch("Lep_E_MaterialGapDown" ,					&Lep_E_MaterialGapDown,					"Lep_E_MaterialGapDown[4]/F"		);
	tree->Branch("Lep_E_MaterialCryoUp" ,					&Lep_E_MaterialCryoUp,					"Lep_E_MaterialCryoUp[4]/F"			);
	tree->Branch("Lep_E_MaterialCryoDown" ,					&Lep_E_MaterialCryoDown,				"Lep_E_MaterialCryoDown[4]/F"		);
	tree->Branch("Lep_E_PileUpUp" ,							&Lep_E_PileUpUp,						"Lep_E_PileUpUp[4]/F"				);
	tree->Branch("Lep_E_PileUpDown" ,						&Lep_E_PileUpDown,						"Lep_E_PileUpDown[4]/F"				);
	tree->Branch("Lep_E_Nom" ,								&Lep_E_Nom,								"Lep_E_Nom[4]/F"					);

	tree->Branch("Lep_E_MaterialAllUp" ,					&Lep_E_MaterialAllUp,					"Lep_E_MaterialAllUp[4]/F"			);
	tree->Branch("Lep_E_MaterialAllDown" ,					&Lep_E_MaterialAllDown,					"Lep_E_MaterialAllDown[4]/F"		);


	tree->Branch("Lep_E_PSUp_Barrel" ,						&Lep_E_PSUp_Barrel,						"Lep_E_PSUp_Barrel[4]/F"					);
	tree->Branch("Lep_E_PSDown_Barrel" ,					&Lep_E_PSDown_Barrel,					"Lep_E_PSDown_Barrel[4]/F"					);
	tree->Branch("Lep_E_S12Up_Barrel" ,						&Lep_E_S12Up_Barrel,					"Lep_E_S12Up_Barrel[4]/F"					);
	tree->Branch("Lep_E_S12Down_Barrel" ,					&Lep_E_S12Down_Barrel,					"Lep_E_S12Down_Barrel[4]/F"					);
	tree->Branch("Lep_E_MatCryoUp_Barrel" ,					&Lep_E_MatCryoUp_Barrel,				"Lep_E_MatCryoUp_Barrel[4]/F"				);
	tree->Branch("Lep_E_MatCryoDown_Barrel" ,				&Lep_E_MatCryoDown_Barrel,				"Lep_E_MatCryoDown_Barrel[4]/F"				);
	tree->Branch("Lep_E_MatCaloUp_Barrel" ,					&Lep_E_MatCaloUp_Barrel,				"Lep_E_MatCaloUp_Barrel[4]/F"				);
	tree->Branch("Lep_E_MatCaloDown_Barrel" ,				&Lep_E_MatCaloDown_Barrel,				"Lep_E_MatCaloDown_Barrel[4]/F"				);
	tree->Branch("Lep_E_LArCalibUp_Barrel" ,				&Lep_E_LArCalibUp_Barrel,				"Lep_E_LArCalibUp_Barrel[4]/F"				);
	tree->Branch("Lep_E_LArCalibDown_Barrel" ,				&Lep_E_LArCalibDown_Barrel,				"Lep_E_LArCalibDown_Barrel[4]/F"			);
	tree->Branch("Lep_E_LArUnconvCalibUp_Barrel" ,			&Lep_E_LArUnconvCalibUp_Barrel,			"Lep_E_LArUnconvCalibUp_Barrel[4]/F"		);
	tree->Branch("Lep_E_LArUnconvCalibDown_Barrel" ,		&Lep_E_LArUnconvCalibDown_Barrel,		"Lep_E_LArUnconvCalibDown_Barrel[4]/F"		);
	tree->Branch("Lep_E_LArElecUnconvUp_Barrel" ,			&Lep_E_LArElecUnconvUp_Barrel,			"Lep_E_LArElecUnconvUp_Barrel[4]/F"			);
	tree->Branch("Lep_E_LArElecUnconvDown_Barrel" ,			&Lep_E_LArElecUnconvDown_Barrel,		"Lep_E_LArElecUnconvDown_Barrel[4]/F"		);

	tree->Branch("Lep_E_PSUp_EC" ,							&Lep_E_PSUp_EC,							"Lep_E_PSUp_EC[4]/F"					);
	tree->Branch("Lep_E_PSDown_EC" ,						&Lep_E_PSDown_EC,						"Lep_E_PSDown_EC[4]/F"					);
	tree->Branch("Lep_E_S12Up_EC" ,							&Lep_E_S12Up_EC,						"Lep_E_S12Up_EC[4]/F"					);
	tree->Branch("Lep_E_S12Down_EC" ,						&Lep_E_S12Down_EC,						"Lep_E_S12Down_EC[4]/F"					);
	tree->Branch("Lep_E_MatCryoUp_EC" ,						&Lep_E_MatCryoUp_EC,					"Lep_E_MatCryoUp_EC[4]/F"				);
	tree->Branch("Lep_E_MatCryoDown_EC" ,					&Lep_E_MatCryoDown_EC,					"Lep_E_MatCryoDown_EC[4]/F"				);
	tree->Branch("Lep_E_MatCaloUp_EC" ,						&Lep_E_MatCaloUp_EC,					"Lep_E_MatCaloUp_EC[4]/F"				);
	tree->Branch("Lep_E_MatCaloDown_EC" ,					&Lep_E_MatCaloDown_EC,					"Lep_E_MatCaloDown_EC[4]/F"				);
	tree->Branch("Lep_E_LArCalibUp_EC" ,					&Lep_E_LArCalibUp_EC,					"Lep_E_LArCalibUp_EC[4]/F"				);
	tree->Branch("Lep_E_LArCalibDown_EC" ,					&Lep_E_LArCalibDown_EC,					"Lep_E_LArCalibDown_EC[4]/F"			);
	tree->Branch("Lep_E_LArUnconvCalibUp_EC" ,				&Lep_E_LArUnconvCalibUp_EC,				"Lep_E_LArUnconvCalibUp_EC[4]/F"		);
	tree->Branch("Lep_E_LArUnconvCalibDown_EC" ,			&Lep_E_LArUnconvCalibDown_EC,			"Lep_E_LArUnconvCalibDown_EC[4]/F"		);
	tree->Branch("Lep_E_LArElecUnconvUp_EC" ,				&Lep_E_LArElecUnconvUp_EC,				"Lep_E_LArElecUnconvUp_EC[4]/F"			);
	tree->Branch("Lep_E_LArElecUnconvDown_EC" ,				&Lep_E_LArElecUnconvDown_EC,			"Lep_E_LArElecUnconvDown_EC[4]/F"		);

	tree->Branch("Lep_E_MatIDUp1" ,							&Lep_E_MatIDUp1,						"Lep_E_MatIDUp1[4]/F"				);
	tree->Branch("Lep_E_MatIDDown1" ,						&Lep_E_MatIDDown1,						"Lep_E_MatIDDown1[4]/F"				);
	tree->Branch("Lep_E_MatIDUp2" ,							&Lep_E_MatIDUp2,						"Lep_E_MatIDUp2[4]/F"				);
	tree->Branch("Lep_E_MatIDDown2" ,						&Lep_E_MatIDDown2,						"Lep_E_MatIDDown2[4]/F"				);
	tree->Branch("Lep_E_MatIDUp3" ,							&Lep_E_MatIDUp3,						"Lep_E_MatIDUp3[4]/F"				);
	tree->Branch("Lep_E_MatIDDown3" ,						&Lep_E_MatIDDown3,						"Lep_E_MatIDDown3[4]/F"				);
	tree->Branch("Lep_E_MatIDUp4" ,							&Lep_E_MatIDUp4,						"Lep_E_MatIDUp4[4]/F"				);
	tree->Branch("Lep_E_MatIDDown4" ,						&Lep_E_MatIDDown4,						"Lep_E_MatIDDown4[4]/F"				);


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
	
	tree->Branch("flagQuad" , 								&tflagQuad ,           					"flagQuad/I"  			);	
}

void OutputTreeSys::fillTree(D3PDReader::Event *event, QuadLepton * higgs, Int_t type, Bool_t isMC)
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

	tm4l_constrained 		= higgs->getMassZMassCons()/1000;
	tm4lerr_constrained 	= higgs->getMassErrZMassCons()/1000;
	tmZ1_constrained 		= higgs->getZ1MassZMassCons()/1000;
	tmZ2_constrained 		= higgs->getZ2MassZMassCons()/1000;

	tweight 		= higgs->weight;
	tweight_corr 	= higgs->weight_corr;
	tweight_lumi 	= higgs->weight_lumi;

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

	tZ1_lepplus_E 		= higgs->getZ1()->getLepPlus()->E_nom/1000;
	tZ1_lepminus_E 		= higgs->getZ1()->getLepNeg()->E_nom/1000;
	tZ2_lepplus_E		= higgs->getZ2()->getLepPlus()->E_nom/1000;
	tZ2_lepminus_E 		= higgs->getZ2()->getLepNeg()->E_nom/1000;

	tZ1_lepplus_id 		= (Int_t) higgs->getZ1()->getLepPlus()->lepID;
	tZ1_lepminus_id 	= (Int_t) higgs->getZ1()->getLepNeg()->lepID;
	tZ2_lepplus_id 		= (Int_t) higgs->getZ2()->getLepPlus()->lepID;
	tZ2_lepminus_id 	= (Int_t) higgs->getZ2()->getLepNeg()->lepID;
	
	tZ1_lepplus_cov_mom  	= higgs->getZ1()->getLepPlus()->covMomErr;
	tZ1_lepminus_cov_mom  	= higgs->getZ1()->getLepNeg()->covMomErr;
	tZ2_lepplus_cov_mom  	= higgs->getZ2()->getLepPlus()->covMomErr;
	tZ2_lepminus_cov_mom  	= higgs->getZ2()->getLepNeg()->covMomErr;
	

	tZ1_lepplus_pt_uncorr_CB  	= higgs->getZ1()->getLepPlus()->cb_pt_unsmeared/1000;
	tZ1_lepminus_pt_uncorr_CB  	= higgs->getZ1()->getLepNeg()->cb_pt_unsmeared/1000;
	tZ2_lepplus_pt_uncorr_CB  	= higgs->getZ2()->getLepPlus()->cb_pt_unsmeared/1000;
	tZ2_lepminus_pt_uncorr_CB  	= higgs->getZ2()->getLepNeg()->cb_pt_unsmeared/1000;
	
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


	vector<ChargedLepton *> currLep = higgs->getLepton();

	for(Int_t i = 0; i < (Int_t) currLep.size(); i++)
	{
		if(currLep[i]->getFlavor() == flavor::Muon)
		{
			Lep_cl_eta[i] = currLep[i]->get4Momentum()->Eta();
		}
		else if(currLep[i]->getFlavor() == flavor::Electron)
		{
			Lep_cl_eta[i] = currLep[i]->GetElectron()->cl_eta();
		}

	}
	// Filling the Sys
	for(Int_t i = 0; i <= doSys::Nom; i++)
	{
		// For each Lepton
		for(Int_t j = 0; j < (Int_t) currLep.size(); j++)
		{
			     if(i == doSys::ZeeStatUp)				Lep_E_ZeeStatUp[j]			= currLep[j]->E_sys[i]/1000;
			else if(i == doSys::ZeeStatDown)			Lep_E_ZeeStatDown[j]		= currLep[j]->E_sys[i]/1000;
			else if(i == doSys::ZeeSystUp)				Lep_E_ZeeSystUp[j]			= currLep[j]->E_sys[i]/1000;
			else if(i == doSys::ZeeSystDown)			Lep_E_ZeeSystDown[j]		= currLep[j]->E_sys[i]/1000;
			else if(i == doSys::ZeeAllUp)				Lep_E_ZeeAllUp[j]			= currLep[j]->E_sys[i]/1000;
			else if(i == doSys::ZeeAllDown)				Lep_E_ZeeAllDown[j]			= currLep[j]->E_sys[i]/1000;
			else if(i == doSys::PSUp)				 	Lep_E_PSUp[j]				= currLep[j]->E_sys[i]/1000;
			else if(i == doSys::PSDown)				 	Lep_E_PSDown[j]				= currLep[j]->E_sys[i]/1000;
			else if(i == doSys::S12Up)				 	Lep_E_S12Up[j]				= currLep[j]->E_sys[i]/1000;
			else if(i == doSys::S12Down)				Lep_E_S12Down[j]			= currLep[j]->E_sys[i]/1000;
			else if(i == doSys::MatIDUp)				Lep_E_MatIDUp[j]			= currLep[j]->E_sys[i]/1000;
			else if(i == doSys::MatIDDown)				Lep_E_MatIDDown[j]			= currLep[j]->E_sys[i]/1000;
			else if(i == doSys::MatCryoUp)				Lep_E_MatCryoUp[j]			= currLep[j]->E_sys[i]/1000;
			else if(i == doSys::MatCryoDown)			Lep_E_MatCryoDown[j]		= currLep[j]->E_sys[i]/1000;
			else if(i == doSys::MatCaloUp)				Lep_E_MatCaloUp[j]			= currLep[j]->E_sys[i]/1000;
			else if(i == doSys::MatCaloDown)			Lep_E_MatCaloDown[j]		= currLep[j]->E_sys[i]/1000;
			else if(i == doSys::LArCalibUp)				Lep_E_LArCalibUp[j]			= currLep[j]->E_sys[i]/1000;
			else if(i == doSys::LArCalibDown)			Lep_E_LArCalibDown[j]		= currLep[j]->E_sys[i]/1000;
			else if(i == doSys::LArUnconvCalibUp)		Lep_E_LArUnconvCalibUp[j]	= currLep[j]->E_sys[i]/1000;
			else if(i == doSys::LArUnconvCalibDown)		Lep_E_LArUnconvCalibDown[j]	= currLep[j]->E_sys[i]/1000;
			else if(i == doSys::LArElecUnconvUp)		Lep_E_LArElecUnconvUp[j]	= currLep[j]->E_sys[i]/1000;
			else if(i == doSys::LArElecUnconvDown)		Lep_E_LArElecUnconvDown[j]	= currLep[j]->E_sys[i]/1000;
			else if(i == doSys::LArElecCalibUp)			Lep_E_LArElecCalibUp[j]		= currLep[j]->E_sys[i]/1000;
			else if(i == doSys::LArElecCalibDown)		Lep_E_LArElecCalibDown[j]	= currLep[j]->E_sys[i]/1000;
			else if(i == doSys::GainUp)					Lep_E_GainUp[j]				= currLep[j]->E_sys[i]/1000;
			else if(i == doSys::GainDown)				Lep_E_GainDown[j]			= currLep[j]->E_sys[i]/1000;
			else if(i == doSys::G4Up)				 	Lep_E_G4Up[j]				= currLep[j]->E_sys[i]/1000;
			else if(i == doSys::G4Down)				 	Lep_E_G4Down[j]				= currLep[j]->E_sys[i]/1000;
			else if(i == doSys::MomentumUp)				Lep_E_MomentumUp[j]			= currLep[j]->E_sys[i]/1000;
			else if(i == doSys::MomentumDown)			Lep_E_MomentumDown[j]		= currLep[j]->E_sys[i]/1000;
			else if(i == doSys::ZSmearingUp)			Lep_E_ZSmearingUp[j]		= currLep[j]->E_sys[i]/1000;
			else if(i == doSys::ZSmearingDown)			Lep_E_ZSmearingDown[j]		= currLep[j]->E_sys[i]/1000;
			else if(i == doSys::SamplingTermUp)			Lep_E_SamplingTermUp[j]		= currLep[j]->E_sys[i]/1000;
			else if(i == doSys::SamplingTermDown)	    Lep_E_SamplingTermDown[j]	= currLep[j]->E_sys[i]/1000;
			else if(i == doSys::MaterialIDUp)			Lep_E_MaterialIDUp[j]		= currLep[j]->E_sys[i]/1000;
			else if(i == doSys::MaterialIDDown)			Lep_E_MaterialIDDown[j]		= currLep[j]->E_sys[i]/1000;
			else if(i == doSys::MaterialCaloUp)			Lep_E_MaterialCaloUp[j]		= currLep[j]->E_sys[i]/1000;
			else if(i == doSys::MaterialCaloDown)	    Lep_E_MaterialCaloDown[j]	= currLep[j]->E_sys[i]/1000;
			else if(i == doSys::MaterialGapUp)			Lep_E_MaterialGapUp[j]		= currLep[j]->E_sys[i]/1000;
			else if(i == doSys::MaterialGapDown)	    Lep_E_MaterialGapDown[j]	= currLep[j]->E_sys[i]/1000;
			else if(i == doSys::MaterialCryoUp)			Lep_E_MaterialCryoUp[j]		= currLep[j]->E_sys[i]/1000;
			else if(i == doSys::MaterialCryoDown)		Lep_E_MaterialCryoDown[j]	= currLep[j]->E_sys[i]/1000;
			else if(i == doSys::PileUpUp)				Lep_E_PileUpUp[j]			= currLep[j]->E_sys[i]/1000;
			else if(i == doSys::PileUpDown)				Lep_E_PileUpDown[j]			= currLep[j]->E_sys[i]/1000;
			else if(i == doSys::Nom)				 	Lep_E_Nom[j]				= currLep[j]->E_sys[i]/1000;

		}
	}

	// For each Lepton
	for(Int_t j = 0; j < (Int_t) currLep.size(); j++)
	{

		// Material ID
		Double_t upRes = 0;
		upRes += pow((Lep_E_MaterialIDUp[j] - Lep_E_Nom[j]),2);
		upRes += pow((Lep_E_MaterialCaloUp[j] - Lep_E_Nom[j]),2);
		upRes += pow((Lep_E_MaterialGapUp[j] - Lep_E_Nom[j]),2);
		upRes += pow((Lep_E_MaterialCryoUp[j] - Lep_E_Nom[j]),2);
		upRes = sqrt(upRes);

		Double_t downRes = 0;
		downRes += pow((Lep_E_MaterialIDDown[j] - Lep_E_Nom[j]),2);
		downRes += pow((Lep_E_MaterialCaloDown[j] - Lep_E_Nom[j]),2);
		downRes += pow((Lep_E_MaterialGapDown[j] - Lep_E_Nom[j]),2);
		downRes += pow((Lep_E_MaterialCryoDown[j] - Lep_E_Nom[j]),2);
		downRes = sqrt(downRes);
		
		Lep_E_MaterialAllUp[j] = Lep_E_Nom[j]  + upRes;
		Lep_E_MaterialAllDown[j] = Lep_E_Nom[j] +  downRes;

		if(fabs(Lep_cl_eta[j]) < 1.4)
		{
			Lep_E_PSUp_Barrel[j]					= Lep_E_PSUp[j];					
			Lep_E_PSDown_Barrel[j]					= Lep_E_PSDown[j];					
			Lep_E_S12Up_Barrel[j]					= Lep_E_S12Up[j];					
			Lep_E_S12Down_Barrel[j]					= Lep_E_S12Down[j];				
			Lep_E_MatCryoUp_Barrel[j]				= Lep_E_MatCryoUp[j];					
			Lep_E_MatCryoDown_Barrel[j]				= Lep_E_MatCryoDown[j];				
			Lep_E_MatCaloUp_Barrel[j]				= Lep_E_MatCaloUp[j];					
			Lep_E_MatCaloDown_Barrel[j]				= Lep_E_MatCaloDown[j];				
			Lep_E_LArUnconvCalibUp_Barrel[j]		= Lep_E_LArUnconvCalibUp[j];					
			Lep_E_LArUnconvCalibDown_Barrel[j]		= Lep_E_LArUnconvCalibDown[j];					
			Lep_E_LArElecUnconvUp_Barrel[j]			= Lep_E_LArElecUnconvUp[j];				
			Lep_E_LArElecUnconvDown_Barrel[j]		= Lep_E_LArElecUnconvDown[j];					
			Lep_E_LArCalibUp_Barrel[j]				= Lep_E_LArCalibUp[j];					
			Lep_E_LArCalibDown_Barrel[j]			= Lep_E_LArCalibDown[j];

			Lep_E_PSUp_EC[j]						= Lep_E_Nom[j];	
			Lep_E_PSDown_EC[j]						= Lep_E_Nom[j];	
			Lep_E_S12Up_EC[j]						= Lep_E_Nom[j];	
			Lep_E_S12Down_EC[j]						= Lep_E_Nom[j];
			Lep_E_MatCryoUp_EC[j]					= Lep_E_Nom[j];		
			Lep_E_MatCryoDown_EC[j]					= Lep_E_Nom[j];	
			Lep_E_MatCaloUp_EC[j]					= Lep_E_Nom[j];		
			Lep_E_MatCaloDown_EC[j]					= Lep_E_Nom[j];	
			Lep_E_LArUnconvCalibUp_EC[j]			= Lep_E_Nom[j];				
			Lep_E_LArUnconvCalibDown_EC[j]			= Lep_E_Nom[j];				
			Lep_E_LArElecUnconvUp_EC[j]				= Lep_E_Nom[j];		
			Lep_E_LArElecUnconvDown_EC[j]			= Lep_E_Nom[j];				
			Lep_E_LArCalibUp_EC[j]					= Lep_E_Nom[j];		
			Lep_E_LArCalibDown_EC[j]				= Lep_E_Nom[j];
		}
		else 
		{
			Lep_E_PSUp_Barrel[j]					= Lep_E_Nom[j];	
			Lep_E_PSDown_Barrel[j]					= Lep_E_Nom[j];	
			Lep_E_S12Up_Barrel[j]					= Lep_E_Nom[j];	
			Lep_E_S12Down_Barrel[j]					= Lep_E_Nom[j];
			Lep_E_MatCryoUp_Barrel[j]				= Lep_E_Nom[j];		
			Lep_E_MatCryoDown_Barrel[j]				= Lep_E_Nom[j];	
			Lep_E_MatCaloUp_Barrel[j]				= Lep_E_Nom[j];		
			Lep_E_MatCaloDown_Barrel[j]				= Lep_E_Nom[j];	
			Lep_E_LArUnconvCalibUp_Barrel[j]		= Lep_E_Nom[j];				
			Lep_E_LArUnconvCalibDown_Barrel[j]		= Lep_E_Nom[j];				
			Lep_E_LArElecUnconvUp_Barrel[j]			= Lep_E_Nom[j];		
			Lep_E_LArElecUnconvDown_Barrel[j]		= Lep_E_Nom[j];				
			Lep_E_LArCalibUp_Barrel[j]				= Lep_E_Nom[j];		
			Lep_E_LArCalibDown_Barrel[j]			= Lep_E_Nom[j];		

			Lep_E_PSUp_EC[j]						= Lep_E_PSUp[j];					
			Lep_E_PSDown_EC[j]						= Lep_E_PSDown[j];					
			Lep_E_S12Up_EC[j]						= Lep_E_S12Up[j];					
			Lep_E_S12Down_EC[j]						= Lep_E_S12Down[j];				
			Lep_E_MatCryoUp_EC[j]					= Lep_E_MatCryoUp[j];					
			Lep_E_MatCryoDown_EC[j]					= Lep_E_MatCryoDown[j];				
			Lep_E_MatCaloUp_EC[j]					= Lep_E_MatCaloUp[j];					
			Lep_E_MatCaloDown_EC[j]					= Lep_E_MatCaloDown[j];				
			Lep_E_LArUnconvCalibUp_EC[j]			= Lep_E_LArUnconvCalibUp[j];					
			Lep_E_LArUnconvCalibDown_EC[j]			= Lep_E_LArUnconvCalibDown[j];					
			Lep_E_LArElecUnconvUp_EC[j]				= Lep_E_LArElecUnconvUp[j];				
			Lep_E_LArElecUnconvDown_EC[j]			= Lep_E_LArElecUnconvDown[j];					
			Lep_E_LArCalibUp_EC[j]					= Lep_E_LArCalibUp[j];					
			Lep_E_LArCalibDown_EC[j]				= Lep_E_LArCalibDown[j];	
		}

		// For material ID
		if(fabs(Lep_cl_eta[j]) < 1.1)
		{
			Lep_E_MatIDUp1[j] = Lep_E_MatIDUp[j];
			Lep_E_MatIDUp2[j] = Lep_E_Nom[j];
			Lep_E_MatIDUp3[j] = Lep_E_Nom[j];
			Lep_E_MatIDUp4[j] = Lep_E_Nom[j];

			Lep_E_MatIDDown1[j] = Lep_E_MatIDDown[j];
			Lep_E_MatIDDown2[j] = Lep_E_Nom[j];
			Lep_E_MatIDDown3[j] = Lep_E_Nom[j];
			Lep_E_MatIDDown4[j] = Lep_E_Nom[j];

			
		}
		else if(fabs(Lep_cl_eta[j]) >= 1.1 && fabs(Lep_cl_eta[j]) < 1.5)
		{
			Lep_E_MatIDUp1[j] = Lep_E_Nom[j];
			Lep_E_MatIDUp2[j] = Lep_E_MatIDUp[j];
			Lep_E_MatIDUp3[j] = Lep_E_Nom[j];
			Lep_E_MatIDUp4[j] = Lep_E_Nom[j];
			
			Lep_E_MatIDDown1[j] = Lep_E_Nom[j];
			Lep_E_MatIDDown2[j] = Lep_E_MatIDDown[j];
			Lep_E_MatIDDown3[j] = Lep_E_Nom[j];
			Lep_E_MatIDDown4[j] = Lep_E_Nom[j];
		}
		else if(fabs(Lep_cl_eta[j]) >= 1.5 && fabs(Lep_cl_eta[j]) < 2.1)
		{
			Lep_E_MatIDUp1[j] = Lep_E_Nom[j];
			Lep_E_MatIDUp2[j] = Lep_E_Nom[j];
			Lep_E_MatIDUp3[j] = Lep_E_MatIDUp[j];
			Lep_E_MatIDUp4[j] = Lep_E_Nom[j];
			
			Lep_E_MatIDDown1[j] = Lep_E_Nom[j];
			Lep_E_MatIDDown2[j] = Lep_E_Nom[j];
			Lep_E_MatIDDown3[j] = Lep_E_MatIDDown[j];
			Lep_E_MatIDDown4[j] = Lep_E_Nom[j];

		}
		else
		{
			Lep_E_MatIDUp1[j] = Lep_E_Nom[j];
			Lep_E_MatIDUp2[j] = Lep_E_Nom[j];
			Lep_E_MatIDUp3[j] = Lep_E_Nom[j];
			Lep_E_MatIDUp4[j] = Lep_E_MatIDUp[j];
			
			Lep_E_MatIDDown1[j] = Lep_E_Nom[j];
			Lep_E_MatIDDown2[j] = Lep_E_Nom[j];
			Lep_E_MatIDDown3[j] = Lep_E_Nom[j];
			Lep_E_MatIDDown4[j] = Lep_E_MatIDDown[j];

		}



	}



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

	// Filling the tree
	if(type == analysisType::Mu4) mu4Tree->Fill();
	else if(type == analysisType::El4) el4Tree->Fill();
	else if(type == analysisType::Mu2El2) mu2el2Tree->Fill();
	else if(type == analysisType::El2Mu2) el2mu2Tree->Fill();
	else {cout<<"OutputTreeSys: fillTree: Analysis type not recognized"<<endl;}
	// Production Category
	if(type == analysisType::Mu4 && higgs->prodChannel == productionChannel::ggF) mu4ggFTree->Fill();
	else if(type == analysisType::El4 && higgs->prodChannel == productionChannel::ggF) el4ggFTree->Fill();
	else if(type == analysisType::Mu2El2 && higgs->prodChannel == productionChannel::ggF) mu2el2ggFTree->Fill();
	else if(type == analysisType::El2Mu2 && higgs->prodChannel == productionChannel::ggF) el2mu2ggFTree->Fill();
	else if(higgs->prodChannel == productionChannel::VBF) VBFTree->Fill();
	else if(higgs->prodChannel == productionChannel::VHLep) VHLepTree->Fill();
	else if(higgs->prodChannel == productionChannel::VHHad) VHHadTree->Fill();	
	else if(higgs->prodChannel == productionChannel::VH) VHTree->Fill();	
	else {cout<<"OutputTreeSys: fillTree: Analysis type not recognized Prod Category"<<endl;}

	// Clearing the vars for the next one
	clearVars();
}

void OutputTreeSys::saveTrees(TString filePath, TH1F* countingHist, TString sampleName)
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


