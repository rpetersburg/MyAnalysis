#include <stdlib.h>
#include <string>
#include "MyAnalysis/ChargedLepton.h"
#include <math.h>

using namespace std;
////////////////////////////////////////////////////////////////////////////////////////
//				Constructor and Destructor
////////////////////////////////////////////////////////////////////////////////////////

ChargedLepton::ChargedLepton(D3PDReader::MuonD3PDObjectElement *tMu, Int_t tType, Double_t tLepEff, Int_t tindex, Bool_t isMC)
{
	mu = tMu;
	type = tType;
	flavor = flavor::Muon;
	smearVal = 1;
	lepEff = tLepEff;
	index = tindex;
	charge = mu->charge();
	elRescale = 0;
	covMomErr = -1;
	lepID = -1;

	// This is done to differentiate between standalone muons and not stand alone muons
	// If this muon is staco, it automatically checks if it is standalone and updates the type
	if(tType == leptonType::MuonStaco)
	{
		if(mu->isStandAloneMuon() == 1)
		{
			type = leptonType::MuonStandAlone;
		}
	}
	
	// Creating the TLorentz Vector
	m_momentum_main = new TLorentzVector();

	Double_t eta = mu->eta();
	Double_t phi = mu->phi();
	Double_t pT = mu->pt();
	Double_t E = mu->E();

	//m_momentum_main->SetPtEtaPhiE(pT, eta, phi, E);
	//m_momentum.SetPtEtaPhiE(pT, eta, phi, E);
	m_momentum_main->SetPtEtaPhiM(pT, eta, phi, pdgMuMass);
	m_momentum.SetPtEtaPhiM(pT, eta, phi, pdgMuMass);
	m_momentumBDT.SetPtEtaPhiM(pT/1000, eta, phi, pdgMuMass/1000);
	// For ms and ID
	m_momentumMS.SetPtEtaPhiM(mu->me_pt, mu->me_eta, mu->me_phi(), pdgMuMass);
	m_momentumID.SetPtEtaPhiM(mu->id_pt, mu->id_eta, mu->id_phi(), pdgMuMass);

	// For Trigger efficiency
	m_momentumTrig = m_momentum;
	// PtCone
	ptCone20 = mu->ptcone20();

	// PtCorrection
	if(type == leptonType::MuonStandAlone) {ptCone20Correction = 0;}
	else {ptCone20Correction = 1./fabs(mu->id_qoverp())*sin(mu->id_theta());}

	// Mass
	lepMass = pdgMuMass;

	// For sys studies
	id_pt_unsmeared = mu->id_pt_unsmeared;  
    me_pt_unsmeared = mu->me_pt_unsmeared;
 	cb_pt_unsmeared = mu->cb_pt_unsmeared;

	for(Int_t i = 0; i < doSys::Nom + 1; i++)
	{
		E_sys[i]= mu->E();
	}
	E_nom = E;

	if(charge > 0) 	lepType = VHLeptonType::muonPlus;
	else 			lepType = VHLeptonType::muonMinus;
	truthParentType = VHTruthType::unknown;	

}
ChargedLepton::ChargedLepton(D3PDReader::ElectronD3PDObjectElement *tEl, Int_t tType, Double_t tLepEff, Double_t tSmearVal, Int_t tindex, Bool_t isMC, Double_t tresolution, Double_t tclPt)
{
	el = tEl;
	type = tType;
	flavor = flavor::Electron;
	lepEff = tLepEff;	
	smearVal = tSmearVal;
	index = tindex;
	charge = el->charge();
	elRescale = 0;
	covMomErr = -1;
	lepID = -1;
	resolutionEl = tresolution;
	clPt = tclPt;
	
	// Creating the TLorentz Vector
	m_momentum_main = new TLorentzVector();

	Double_t phi_trk = el->trackphi();
	Double_t eta_trk = el->tracketa();
	Double_t E = el->cl_E();
	Double_t pT = E/cosh(eta_trk);
	
	m_momentum_main->SetPtEtaPhiM(pT, eta_trk, phi_trk, pdgElMass);
	m_momentum.SetPtEtaPhiM(pT, eta_trk, phi_trk, pdgElMass);
	m_momentumBDT.SetPtEtaPhiM(pT/1000, eta_trk, phi_trk, pdgElMass/1000);
	m_momentumMS = m_momentum;
	m_momentumID = m_momentum;
	//m_momentum_main->SetPtEtaPhiE(pT, eta_trk, phi_trk, E);
	//m_momentum.SetPtEtaPhiE(pT, eta_trk, phi_trk, E);
	
	// TLorentz Vector for Trigger Efficiency
	Double_t phi_cl = el->cl_phi();
	Double_t eta_cl = el->cl_eta();
	//Double_t E_cl = el->cl_E();
	//Double_t pT_cl = E/cosh(eta_cl);
	//m_momentumTrig.SetPtEtaPhiE(pT, eta_cl, phi_cl, E);
	m_momentumTrig.SetPtEtaPhiM(el->bfEP_pT, eta_cl, phi_cl, pdgElMass);

	// PtCone
	ptCone20 = el->ptcone20();

	ptCone20Correction = el->trackpt();

	// Mass
	lepMass = pdgElMass;
	
	
	// For Sys
	id_pt_unsmeared = el->pT_unsmeared;  
    me_pt_unsmeared = el->pT_unsmeared;
    cb_pt_unsmeared = el->pT_unsmeared;

	for(Int_t i = 0; i < doSys::Nom + 1; i++)
	{
		E_sys[i]= el->E_sysVar[i];
	}
	E_nom = E;
	
	if(charge > 0) 	lepType = VHLeptonType::electronPlus;
	else 			lepType = VHLeptonType::electronMinus;
	truthParentType = VHTruthType::unknown;	

}

ChargedLepton::ChargedLepton(D3PDReader::JetD3PDObjectElement *tJets, Int_t tType, Int_t tindex, Bool_t isMC)
{
	jets = tJets;
	type = tType;
	index = tindex;
	flavor = flavor::Jet;
	charge = 0;
	elRescale = 0;
	covMomErr = -1;
	
	// Creating the TLorentz Vector
	m_momentum_main = new TLorentzVector();

//	Double_t eta = jets->emscale_eta();
	Double_t eta = jets->eta();	
	Double_t phi = jets->phi();
	Double_t pT = jets->pt();
//	Double_t E = jets->E();
	Double_t M = jets->m();
//	m_momentum_main->SetPtEtaPhiE(pT, eta, phi, E);
//	m_momentum.SetPtEtaPhiE(pT, eta, phi, E);
	// To compare with valerio
	m_momentum_main->SetPtEtaPhiM(pT, eta, phi, M);
	m_momentum.SetPtEtaPhiM(pT, eta, phi, M);
	m_momentumBDT.SetPtEtaPhiM(pT/1000, eta, phi, M/1000);

<<<<<<< HEAD
=======

	//cout<<"Jet M: "<<M<<endl;
	//cout<<"TVL JET M: "<<m_momentum_main->M()<<endl;
>>>>>>> bc7b9ddaf72f0a41dfe1bb5d9068cc4b03444c0d
	// PtCone
	ptCone20 = -1;

	// Mass
	lepMass = -1;

	lepType = VHLeptonType::unknown;
	truthParentType = VHTruthType::unknown;	

}

ChargedLepton::~ChargedLepton()
{
	if(m_momentum_main != 0) delete m_momentum_main;
}

void ChargedLepton::FillCovMatrix(Int_t runNumber_sf)
{
	if(elRescale == 0)
	{
		cout<<"ChargedLepton: ERROR = Energy rescaler not set"<<endl;
		return;
	}

	if(flavor == flavor::Electron)
	{
		Double_t E = el->cl_E();
		Double_t Eta = el->cl_eta();
		Double_t res_electron = elRescale->resolution(E, Eta,PATCore::ParticleType::Electron, true)*E;

		// Overwrite it with the outside resolution if given
		if(resolutionEl != -1 ) res_electron = resolutionEl;
		// Momentum Error Term
		covMomErr =	pow((res_electron),2) /(1000*1000); 
		//covMomErr = res_electron/1000;
		
		TMatrixD tmp1 = ZMassConstraint::getCovarianceTMatrixDd0z0PhiThetaPElectron(res_electron,
	                                                                              el->trackcov_d0(),  
	                                                                              el->trackcov_z0(),  
	                                                                              el->trackcov_phi(),  
	                                                                              el->trackcov_theta(),  
	                                                           	                  el->trackcov_d0_z0(),  
	                                                                              el->trackcov_d0_phi(),  
	                                                                              el->trackcov_d0_theta(),
	                                                                              el->trackcov_z0_phi(),
	                                                                              el->trackcov_z0_theta(),  
	                                                                              el->trackcov_phi_theta());
		errorMatrix.ResizeTo(tmp1);
		errorMatrix = tmp1;

		errorMatrixHep = ZMassConstraint::getCovarianceMatrixd0z0PhiThetaPElectron(res_electron,
	                                                                              el->trackcov_d0(),  
	                                                                              el->trackcov_z0(),  
	                                                                              el->trackcov_phi(),  
	                                                                              el->trackcov_theta(),  
	                                                           	                  el->trackcov_d0_z0(),  
	                                                                              el->trackcov_d0_phi(),  
	                                                                              el->trackcov_d0_theta(),
	                                                                              el->trackcov_z0_phi(),
	                                                                              el->trackcov_z0_theta(),  
	                                                                              el->trackcov_phi_theta());


		//cout<<"el_E_Rol"<<res_electron<<endl;
		//cout<<"trackcov_d0 "<<  el->trackcov_d0()<<endl;
		//cout<<"trackcov_z0 "<<  el->trackcov_z0() <<endl;
		//cout<<"trackcov_phi "<< el->trackcov_phi() <<endl;
		//cout<<"trackcov_theta "<< el->trackcov_theta() <<endl;
		//cout<<"trackcov_d0_z0 "<< el->trackcov_d0_z0() <<endl;
		//cout<<"trackcov_d0_phi "<< el->trackcov_d0_phi() <<endl;
		//cout<<"trackcov_d0_theta "<< el->trackcov_d0_theta() <<endl;
		//cout<<"trackcov_z0_phi "<< el->trackcov_z0_phi() <<endl;
		//cout<<"trackcov_z0_theta "<< el->trackcov_z0_theta() <<endl;
		//cout<<"trackcov_phi_theta "<< el->trackcov_phi_theta() <<endl;
		errorMatrixMS.ResizeTo(tmp1);
		errorMatrixMS = errorMatrix;
		errorMatrixHepMS = errorMatrixHep;
		
		errorMatrixID.ResizeTo(tmp1);
		errorMatrixID = errorMatrix;
		errorMatrixHepID = errorMatrixHep;


	}
	if(flavor == flavor::Muon)
	{
		TMatrixD tmp1 =ZMassConstraint::getCovarianceTMatrixDd0z0PhiThetaPMuon(m_momentum.P(), 
                                                                             mu->cov_d0_exPV(),  
                                                                             mu->cov_z0_exPV(),  
                                                                             mu->cov_phi_exPV(),  
                                                                             mu->cov_theta_exPV(),  
                                                                             mu->cov_qoverp_exPV(), 
                                                                             mu->cov_d0_z0_exPV(),  
                                                                             mu->cov_d0_phi_exPV(),  
                                                                             mu->cov_d0_theta_exPV(),  
                                                                             mu->cov_d0_qoverp_exPV(),  
                                                                             mu->cov_z0_phi_exPV(),
                                                                             mu->cov_z0_theta_exPV(),  
                                                                             mu->cov_z0_qoverp_exPV(),  
                                                                             mu->cov_phi_theta_exPV(),  
                                                                             mu->cov_phi_qoverp_exPV(), 
                                                                             mu->cov_theta_qoverp_exPV());
		errorMatrix.ResizeTo(tmp1);
		errorMatrix = tmp1;

		// Momentum Error Term
		covMomErr = (mu->cov_qoverp_exPV())*1000*1000;
		
		errorMatrixHep = ZMassConstraint::getCovarianceMatrixd0z0PhiThetaPMuon(m_momentum.P(), 
                                                                             mu->cov_d0_exPV(),  
                                                                             mu->cov_z0_exPV(),  
                                                                             mu->cov_phi_exPV(),  
                                                                             mu->cov_theta_exPV(),  
                                                                             mu->cov_qoverp_exPV(), 
                                                                             mu->cov_d0_z0_exPV(),  
                                                                             mu->cov_d0_phi_exPV(),  
                                                                             mu->cov_d0_theta_exPV(),  
                                                                             mu->cov_d0_qoverp_exPV(),  
                                                                             mu->cov_z0_phi_exPV(),
                                                                             mu->cov_z0_theta_exPV(),  
                                                                             mu->cov_z0_qoverp_exPV(),  
                                                                             mu->cov_phi_theta_exPV(),  
                                                                             mu->cov_phi_qoverp_exPV(), 
                                                                             mu->cov_theta_qoverp_exPV());
		
		
		//------------------------
		// For MS
		TMatrixD tmp1MS =ZMassConstraint::getCovarianceTMatrixDd0z0PhiThetaPMuon(m_momentumMS.P(), 
                                                                             	mu->me_cov_d0_exPV(),  
                                                                             	mu->me_cov_z0_exPV(),  
                                                                             	mu->me_cov_phi_exPV(),  
                                                                             	mu->me_cov_theta_exPV(),  
                                                                             	mu->me_cov_qoverp_exPV(), 
                                                                             	mu->me_cov_d0_z0_exPV(),  
                                                                             	mu->me_cov_d0_phi_exPV(),  
                                                                             	mu->me_cov_d0_theta_exPV(),  
                                                                             	mu->me_cov_d0_qoverp_exPV(),  
                                                                             	mu->me_cov_z0_phi_exPV(),
                                                                             	mu->me_cov_z0_theta_exPV(),  
                                                                             	mu->me_cov_z0_qoverp_exPV(),  
                                                                             	mu->me_cov_phi_theta_exPV(),  
                                                                             	mu->me_cov_phi_qoverp_exPV(), 
                                                                             	mu->me_cov_theta_qoverp_exPV());
		errorMatrixMS.ResizeTo(tmp1MS);
		errorMatrixMS = tmp1MS;

		
		errorMatrixHepMS = ZMassConstraint::getCovarianceMatrixd0z0PhiThetaPMuon(m_momentumMS.P(), 
                                                                             mu->me_cov_d0_exPV(),  
                                                                             mu->me_cov_z0_exPV(),  
                                                                             mu->me_cov_phi_exPV(),  
                                                                             mu->me_cov_theta_exPV(),  
                                                                             mu->me_cov_qoverp_exPV(), 
                                                                             mu->me_cov_d0_z0_exPV(),  
                                                                             mu->me_cov_d0_phi_exPV(),  
                                                                             mu->me_cov_d0_theta_exPV(),  
                                                                             mu->me_cov_d0_qoverp_exPV(),  
                                                                             mu->me_cov_z0_phi_exPV(),
                                                                             mu->me_cov_z0_theta_exPV(),  
                                                                             mu->me_cov_z0_qoverp_exPV(),  
                                                                             mu->me_cov_phi_theta_exPV(),  
                                                                             mu->me_cov_phi_qoverp_exPV(), 
                                                                             mu->me_cov_theta_qoverp_exPV());


		//------------------------
		// For ID 
		TMatrixD tmp1ID =ZMassConstraint::getCovarianceTMatrixDd0z0PhiThetaPMuon(m_momentumID.P(), 
                                                                             	mu->id_cov_d0_exPV(),  
                                                                             	mu->id_cov_z0_exPV(),  
                                                                             	mu->id_cov_phi_exPV(),  
                                                                             	mu->id_cov_theta_exPV(),  
                                                                             	mu->id_cov_qoverp_exPV(), 
                                                                             	mu->id_cov_d0_z0_exPV(),  
                                                                             	mu->id_cov_d0_phi_exPV(),  
                                                                             	mu->id_cov_d0_theta_exPV(),  
                                                                             	mu->id_cov_d0_qoverp_exPV(),  
                                                                             	mu->id_cov_z0_phi_exPV(),
                                                                             	mu->id_cov_z0_theta_exPV(),  
                                                                             	mu->id_cov_z0_qoverp_exPV(),  
                                                                             	mu->id_cov_phi_theta_exPV(),  
                                                                             	mu->id_cov_phi_qoverp_exPV(), 
                                                                             	mu->id_cov_theta_qoverp_exPV());
		errorMatrixID.ResizeTo(tmp1ID);
		errorMatrixID = tmp1ID;

		
		errorMatrixHepID = ZMassConstraint::getCovarianceMatrixd0z0PhiThetaPMuon(m_momentumID.P(), 
                                                                             mu->id_cov_d0_exPV(),  
                                                                             mu->id_cov_z0_exPV(),  
                                                                             mu->id_cov_phi_exPV(),  
                                                                             mu->id_cov_theta_exPV(),  
                                                                             mu->id_cov_qoverp_exPV(), 
                                                                             mu->id_cov_d0_z0_exPV(),  
                                                                             mu->id_cov_d0_phi_exPV(),  
                                                                             mu->id_cov_d0_theta_exPV(),  
                                                                             mu->id_cov_d0_qoverp_exPV(),  
                                                                             mu->id_cov_z0_phi_exPV(),
                                                                             mu->id_cov_z0_theta_exPV(),  
                                                                             mu->id_cov_z0_qoverp_exPV(),  
                                                                             mu->id_cov_phi_theta_exPV(),  
                                                                             mu->id_cov_phi_qoverp_exPV(), 
                                                                             mu->id_cov_theta_qoverp_exPV());

		if(mu->isStandAloneMuon())
		{
			errorMatrixHepID = errorMatrixHep;
			errorMatrixID = errorMatrix;
		}
		if(mu->isSegmentTaggedMuon() || mu->isCaloMuonId())
		{
			errorMatrixHepMS = errorMatrixHep;
			errorMatrixMS = errorMatrix;
		}

<<<<<<< HEAD
		//cout<<std::setprecision(10)<<std::scientific<<"CB"<<endl;
=======
		//cout<<std::setprecision(10)<<std::scientific<<endl<<"CB"<<endl;
>>>>>>> bc7b9ddaf72f0a41dfe1bb5d9068cc4b03444c0d
		//cout<<std::setprecision(10)<<std::fixed<<"mu_p "<<m_momentum.P()<<endl;
		//cout<<std::setprecision(10)<<std::scientific<<"cov_d0_exPV "<<  mu->cov_d0_exPV()<<endl;
		//cout<<std::setprecision(10)<<std::scientific<<"cov_z0_exPV "<< mu->cov_z0_exPV() <<endl;
		//cout<<std::setprecision(10)<<std::scientific<<"cov_phi_exPV "<< mu->cov_phi_exPV() <<endl;
		//cout<<std::setprecision(10)<<std::scientific<<"cov_theta_exPV "<< mu->cov_theta_exPV() <<endl;
		//cout<<std::setprecision(10)<<std::scientific<<"cov_qoverp_exPV "<< mu->cov_qoverp_exPV() <<endl;
		//cout<<std::setprecision(10)<<std::scientific<<"cov_d0_z0_exPV "<< mu->cov_d0_z0_exPV() <<endl;
		//cout<<std::setprecision(10)<<std::scientific<<"cov_d0_phi_exPV "<< mu->cov_d0_phi_exPV() <<endl;
		//cout<<std::setprecision(10)<<std::scientific<<"cov_d0_theta_exPV"<<mu->cov_d0_theta_exPV()<<endl;
		//cout<<std::setprecision(10)<<std::scientific<<"cov_d0_qoverp_exPV "<< mu->cov_d0_qoverp_exPV() <<endl;
		//cout<<std::setprecision(10)<<std::scientific<<"cov_z0_phi_exPV "<< mu->cov_z0_phi_exPV() <<endl;
		//cout<<std::setprecision(10)<<std::scientific<<"cov_z0_theta_exPV "<<mu->cov_z0_theta_exPV()<<endl;
		//cout<<std::setprecision(10)<<std::scientific<<"cov_z0_qoverp_exPV "<< mu->cov_z0_qoverp_exPV() <<endl;
		//cout<<std::setprecision(10)<<std::scientific<<"cov_phi_theta_exPV "<< mu->cov_phi_theta_exPV() <<endl;
		//cout<<std::setprecision(10)<<std::scientific<<"cov_phi_qoverp_exPV "<< mu->cov_phi_qoverp_exPV() <<endl;
		//cout<<std::setprecision(10)<<std::scientific<<"cov_theta_qoverp_exPV "<< mu->cov_theta_qoverp_exPV() <<endl;

		//cout<<std::setprecision(10)<<std::scientific<<endl<<"MS"<<endl;
		//cout<<std::setprecision(10)<<std::fixed<<"mu_p "<<m_momentumMS.P()<<endl;
		//cout<<std::setprecision(10)<<std::scientific<<"ms_cov_d0_exPV "<<  mu->me_cov_d0_exPV()<<endl;
		//cout<<std::setprecision(10)<<std::scientific<<"ms_cov_z0_exPV "<< mu->me_cov_z0_exPV() <<endl;
		//cout<<std::setprecision(10)<<std::scientific<<"ms_cov_phi_exPV "<< mu->me_cov_phi_exPV() <<endl;
		//cout<<std::setprecision(10)<<std::scientific<<"ms_cov_theta_exPV "<< mu->me_cov_theta_exPV() <<endl;
		//cout<<std::setprecision(10)<<std::scientific<<"ms_cov_qoverp_exPV "<< mu->me_cov_qoverp_exPV() <<endl;
		//cout<<std::setprecision(10)<<std::scientific<<"ms_cov_d0_z0_exPV "<< mu->me_cov_d0_z0_exPV() <<endl;
		//cout<<std::setprecision(10)<<std::scientific<<"ms_cov_d0_phi_exPV "<< mu->me_cov_d0_phi_exPV() <<endl;
		//cout<<std::setprecision(10)<<std::scientific<<"ms_cov_d0_theta_exPV"<<mu->me_cov_d0_theta_exPV()<<endl;
		//cout<<std::setprecision(10)<<std::scientific<<"ms_cov_d0_qoverp_exPV "<< mu->me_cov_d0_qoverp_exPV() <<endl;
		//cout<<std::setprecision(10)<<std::scientific<<"ms_cov_z0_phi_exPV "<< mu->me_cov_z0_phi_exPV() <<endl;
		//cout<<std::setprecision(10)<<std::scientific<<"ms_cov_z0_theta_exPV "<<mu->me_cov_z0_theta_exPV()<<endl;
		//cout<<std::setprecision(10)<<std::scientific<<"ms_cov_z0_qoverp_exPV "<< mu->me_cov_z0_qoverp_exPV() <<endl;
		//cout<<std::setprecision(10)<<std::scientific<<"ms_cov_phi_theta_exPV "<< mu->me_cov_phi_theta_exPV() <<endl;
		//cout<<std::setprecision(10)<<std::scientific<<"ms_cov_phi_qoverp_exPV "<< mu->me_cov_phi_qoverp_exPV() <<endl;
		//cout<<std::setprecision(10)<<std::scientific<<"ms_cov_theta_qoverp_exPV "<< mu->me_cov_theta_qoverp_exPV() <<endl;


		//cout<<std::setprecision(10)<<std::scientific<<endl<<"ID"<<endl;
		//cout<<std::setprecision(10)<<std::fixed<<"mu_p "<<m_momentumID.P()<<endl;
		//cout<<std::setprecision(10)<<std::scientific<<"id_cov_d0_exPV "<<  mu->id_cov_d0_exPV()<<endl;
		//cout<<std::setprecision(10)<<std::scientific<<"id_cov_z0_exPV "<< mu->id_cov_z0_exPV() <<endl;
		//cout<<std::setprecision(10)<<std::scientific<<"id_cov_phi_exPV "<< mu->id_cov_phi_exPV() <<endl;
		//cout<<std::setprecision(10)<<std::scientific<<"id_cov_theta_exPV "<< mu->id_cov_theta_exPV() <<endl;
		//cout<<std::setprecision(10)<<std::scientific<<"id_cov_qoverp_exPV "<< mu->id_cov_qoverp_exPV() <<endl;
		//cout<<std::setprecision(10)<<std::scientific<<"id_cov_d0_z0_exPV "<< mu->id_cov_d0_z0_exPV() <<endl;
		//cout<<std::setprecision(10)<<std::scientific<<"id_cov_d0_phi_exPV "<< mu->id_cov_d0_phi_exPV() <<endl;
		//cout<<std::setprecision(10)<<std::scientific<<"id_cov_d0_theta_exPV"<<mu->id_cov_d0_theta_exPV()<<endl;
		//cout<<std::setprecision(10)<<std::scientific<<"id_cov_d0_qoverp_exPV "<< mu->id_cov_d0_qoverp_exPV() <<endl;
		//cout<<std::setprecision(10)<<std::scientific<<"id_cov_z0_phi_exPV "<< mu->id_cov_z0_phi_exPV() <<endl;
		//cout<<std::setprecision(10)<<std::scientific<<"id_cov_z0_theta_exPV "<<mu->id_cov_z0_theta_exPV()<<endl;
		//cout<<std::setprecision(10)<<std::scientific<<"id_cov_z0_qoverp_exPV "<< mu->id_cov_z0_qoverp_exPV() <<endl;
		//cout<<std::setprecision(10)<<std::scientific<<"id_cov_phi_theta_exPV "<< mu->id_cov_phi_theta_exPV() <<endl;
		//cout<<std::setprecision(10)<<std::scientific<<"id_cov_phi_qoverp_exPV "<< mu->id_cov_phi_qoverp_exPV() <<endl;
		//cout<<std::setprecision(10)<<std::scientific<<"id_cov_theta_qoverp_exPV "<< mu->id_cov_theta_qoverp_exPV() <<endl;

	}
}

