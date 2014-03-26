#include <stdlib.h>
#include <string>
#include "MyAnalysis/QuadLepton.h"
#include <iostream>


using namespace std;
////////////////////////////////////////////////////////////////////////////////////////
//							Constructor and Destructor
////////////////////////////////////////////////////////////////////////////////////////

QuadLepton::QuadLepton()
{
	Reset();
}

QuadLepton::QuadLepton(DiLepton *ZOne, DiLepton *ZTwo)
{
	Reset();
	Set(ZOne, ZTwo);
}
QuadLepton::~QuadLepton()
{
	if(m_momentum_main != 0) delete m_momentum_main;
	
	leptonInfo.clear();
	
	leptonLorentz.clear();
	leptonCovMatrix.clear();
	leptonHepCovMatrix.clear();	
	
	leptonFSRLorentz.clear();
	leptonFSRCovMatrix.clear();
	leptonHepFSRCovMatrix.clear();
	
	leptonZMassFSRLorentz.clear();
	leptonZMassCovMatrix.clear();

	leptonLorentzID.clear();
	leptonCovMatrixID.clear();
	leptonHepCovMatrixID.clear();
	leptonLorentzMS.clear();
	leptonCovMatrixMS.clear();
	leptonHepCovMatrixMS.clear();

	leptonFSRLorentzID.clear();
	leptonFSRCovMatrixID.clear();
	leptonHepFSRCovMatrixID.clear();
	leptonFSRLorentzMS.clear();
	leptonFSRCovMatrixMS.clear();
	leptonHepFSRCovMatrixMS.clear();

	leptonZMassFSRLorentzID.clear();
	leptonZMassCovMatrixID.clear();
	leptonZMassFSRLorentzMS.clear();
	leptonZMassCovMatrixMS.clear();
	
	looseElectron.clear();
	trackIso.clear();
	caloIso.clear();
	d0Sig.clear();
	valTrackIso.clear();
	valCaloIso.clear();
	valD0Sig.clear();


}

////////////////////////////////////////////////////////////////////////////////////////
//									Helpers				
////////////////////////////////////////////////////////////////////////////////////////
// Resets all the varibles
void QuadLepton::Reset()
{
	quadType = -999; 
	neutral = -999; 

	m_momentum_main = new TLorentzVector();

	Z1 = 0;
	Z2 = 0;

	// Just in case
	mass = -999;
	massErr = -999;
	massID = -999;
	massErrID = -999;	
	massMS = -999;
	massErrMS = -999;
 	massFSR = -999;
 	massErrFSR = -999;
	fsrType = -999;
	massZMassCons = -999;
	massErrZmassCons = -999;
	massZMassConsMS = -999;
	massErrZmassConsMS = -999;
	massZMassConsID = -999;
	massErrZmassConsID = -999;
	massZ1ZMassCons = -999;
	massZ2ZMassCons = -999;
	prodChannel = -999;
	cthstr = -999;
	phi1 = -999;
	cth1 = -999;
	cth2 = -999;
	phi = -999;
	weight = -999;
	weight_corr = -999;
	weight_lumi = -999;
	zVertWeight = -999;
	pileupWeight = -999;
	ttBarVeto = -999;
	ggFWeight = -999;
	JHUWeight = -999;
	zzOverlapVeto = -999;
	zBBVeto = -999;
	trigEff = -999;
	lepEff = -999;
	mcEventWeight = -999;
	crossSection = -999;
	branchRatio = -999;
	lumi = -999;
	flagQuad = -999;

	m4l_truth_born = -999; 
	mZ1_truth_born = -999;
	mZ2_truth_born = -999;
	m4l_truth_bare = -999; 
	mZ1_truth_bare = -999; 
	mZ2_truth_bare = -999;
	m4l_truth_dressed = -999; 
	mZ1_truth_dressed = -999; 
	mZ2_truth_dressed = -999;
	
	leadingJet = 0;
	subLeadingJet = 0;
	thirdJet = 0;
	n_jets = -999;
	dijet_invmass = -999;
	dijet_deltaeta = -999;
	leading_jet_pt = -999;
	leading_jet_eta = -999;
	subleading_jet_pt = -999;
	n_jets_truth_bare = -999;
	leading_jet_pt_truth_bare = -999;

	BDT_discriminant_VBF = -999;
	BDT_discriminant_HadVH = -999;

	KD_discriminant = -999;
	BDT_discriminant = -999;
	BDTGuass_discriminant = -999;
	ptSysupFac = -999;
	ptSysdownFac = -999;

	leadingExtraLep = 0;
	subleadingExtraLep = 0;

	npv = -999;

	calib = -999;

	n_jets_fid = -999;
	leading_jet_pt_fid = -999;
	n_jets_truth_fid = -999;
	leading_jet_pt_truth_fid = -999;
	BCHCutMedium = -999;
	BCHCutTight	 = -999;


}

// Sets all the variable and order comes from outside
void QuadLepton::Set(DiLepton *ZOne, DiLepton *ZTwo)
{
	if(ZOne->getFlavor() == flavor::Muon)
	{
		if(ZTwo->getFlavor() == flavor::Muon)
			quadType = quadType::Mu4;

		else 
			quadType = quadType::Mu2El2;
	}
	else if(ZOne->getFlavor() == flavor::Electron)
	{
		if(ZTwo->getFlavor() == flavor::Electron)
			quadType = quadType::El4;

		else 
			quadType = quadType::El2Mu2;
	}

	quadTypeBR = quadType;

	neutral = (ZOne->IsNeutral() & ZTwo->IsNeutral());

	*m_momentum_main = *ZOne->get4Momentum() + *ZTwo->get4Momentum();
	m_momentum = ZOne->get4MomentumNoP() + ZTwo->get4MomentumNoP();

	Z1 = ZOne;
	Z2 = ZTwo;

	// Filling the vector for the leptons
	// Order is important
	leptonInfo.push_back(ZOne->getLepPlus());	
	leptonInfo.push_back(ZOne->getLepNeg());	
	leptonInfo.push_back(ZTwo->getLepPlus());	
	leptonInfo.push_back(ZTwo->getLepNeg());

	leptonLorentz.push_back(ZOne->getLepPlus()->get4MomentumNoP());
	leptonLorentz.push_back(ZOne->getLepNeg()->get4MomentumNoP());
	leptonLorentz.push_back(ZTwo->getLepPlus()->get4MomentumNoP());
	leptonLorentz.push_back(ZTwo->getLepNeg()->get4MomentumNoP());

	leptonLorentzMS.push_back(ZOne->getLepPlus()->get4MomentumNoP(muonType::MS));
	leptonLorentzMS.push_back(ZOne->getLepNeg()->get4MomentumNoP(muonType::MS));
	leptonLorentzMS.push_back(ZTwo->getLepPlus()->get4MomentumNoP(muonType::MS));
	leptonLorentzMS.push_back(ZTwo->getLepNeg()->get4MomentumNoP(muonType::MS));

	leptonLorentzID.push_back(ZOne->getLepPlus()->get4MomentumNoP(muonType::ID));
	leptonLorentzID.push_back(ZOne->getLepNeg()->get4MomentumNoP(muonType::ID));
	leptonLorentzID.push_back(ZTwo->getLepPlus()->get4MomentumNoP(muonType::ID));
	leptonLorentzID.push_back(ZTwo->getLepNeg()->get4MomentumNoP(muonType::ID));


	leptonFSRLorentz=leptonLorentz;
	leptonFSRLorentzID=leptonLorentzID;
	leptonFSRLorentzMS=leptonLorentzMS;
}

void QuadLepton::SetElRescale(AtlasRoot::egammaEnergyCorrectionTool *telRescale)
{
	elRescale = telRescale;
	Z1->getLepPlus()->SetElRescale(telRescale);
	Z1->getLepNeg()->SetElRescale(telRescale);
	Z2->getLepPlus()->SetElRescale(telRescale);
	Z2->getLepNeg()->SetElRescale(telRescale);
}

void QuadLepton::FillCovMatrix(Int_t runNumber_sf)
{
	Z1->getLepPlus()->FillCovMatrix(runNumber_sf);
	Z1->getLepNeg()->FillCovMatrix(runNumber_sf);
	Z2->getLepPlus()->FillCovMatrix(runNumber_sf);
	Z2->getLepNeg()->FillCovMatrix(runNumber_sf);
	
	// Order is important
	leptonCovMatrix.push_back(Z1->getLepPlus()->getCovMatrix());
	leptonCovMatrix.push_back(Z1->getLepNeg()->getCovMatrix());
	leptonCovMatrix.push_back(Z2->getLepPlus()->getCovMatrix());
	leptonCovMatrix.push_back(Z2->getLepNeg()->getCovMatrix());

	leptonHepCovMatrix.push_back(Z1->getLepPlus()->getHepCovMatrix());
	leptonHepCovMatrix.push_back(Z1->getLepNeg()->getHepCovMatrix());
	leptonHepCovMatrix.push_back(Z2->getLepPlus()->getHepCovMatrix());
	leptonHepCovMatrix.push_back(Z2->getLepNeg()->getHepCovMatrix());

	leptonFSRCovMatrix = leptonCovMatrix;
	leptonHepFSRCovMatrix = leptonHepCovMatrix;

	// For MS
	leptonCovMatrixMS.push_back(Z1->getLepPlus()->getCovMatrix(muonType::MS));
	leptonCovMatrixMS.push_back(Z1->getLepNeg()->getCovMatrix(muonType::MS));
	leptonCovMatrixMS.push_back(Z2->getLepPlus()->getCovMatrix(muonType::MS));
	leptonCovMatrixMS.push_back(Z2->getLepNeg()->getCovMatrix(muonType::MS));

	leptonHepCovMatrixMS.push_back(Z1->getLepPlus()->getHepCovMatrix(muonType::MS));
	leptonHepCovMatrixMS.push_back(Z1->getLepNeg()->getHepCovMatrix(muonType::MS));
	leptonHepCovMatrixMS.push_back(Z2->getLepPlus()->getHepCovMatrix(muonType::MS));
	leptonHepCovMatrixMS.push_back(Z2->getLepNeg()->getHepCovMatrix(muonType::MS));

	leptonFSRCovMatrixMS = leptonCovMatrixMS;
	leptonHepFSRCovMatrixMS = leptonHepCovMatrixMS;


	// For ID
	leptonCovMatrixID.push_back(Z1->getLepPlus()->getCovMatrix(muonType::ID));
	leptonCovMatrixID.push_back(Z1->getLepNeg()->getCovMatrix(muonType::ID));
	leptonCovMatrixID.push_back(Z2->getLepPlus()->getCovMatrix(muonType::ID));
	leptonCovMatrixID.push_back(Z2->getLepNeg()->getCovMatrix(muonType::ID));

	leptonHepCovMatrixID.push_back(Z1->getLepPlus()->getHepCovMatrix(muonType::ID));
	leptonHepCovMatrixID.push_back(Z1->getLepNeg()->getHepCovMatrix(muonType::ID));
	leptonHepCovMatrixID.push_back(Z2->getLepPlus()->getHepCovMatrix(muonType::ID));
	leptonHepCovMatrixID.push_back(Z2->getLepNeg()->getHepCovMatrix(muonType::ID));

	leptonFSRCovMatrixID = leptonCovMatrixID;
	leptonHepFSRCovMatrixID = leptonHepCovMatrixID;

}

void QuadLepton::FillFSRCorr(TLorentzVector momFSR, Int_t runNumber_sf, Bool_t isZ1, Bool_t isZ2, PATCore::ParticleType::Type type)
{
	fsr_momentum = momFSR;
	leptonFSRLorentz.push_back(momFSR);
	leptonFSRLorentzID.push_back(momFSR);
	leptonFSRLorentzMS.push_back(momFSR);

	double energyResolution = elRescale->resolution(momFSR.E(),momFSR.Eta(),type, true)*momFSR.E(); 

    TMatrixD photonCovarianceTMatrixD=ZMassConstraint::getCovarianceTMatrixDd0z0PhiThetaPElectron(energyResolution,           
                                                                                                  0.000001,  
                                                                                                  0.000001,  
                                                                                                  0.000001,  
                                                                                                  0.000001,  
                                                                                                  0.0,  
                                                                                                  0.0,  
                                                                                                  0.0,
                                                                                                  0.0,
                                                                                                  0.0,  
                                                                                                  0.0);

	HepMatrix photonCovarianceHepTMatrixD=ZMassConstraint::getCovarianceMatrixd0z0PhiThetaPElectron(energyResolution,           
                                                                                                  0.000001,  
                                                                                                  0.000001,  
                                                                                                  0.000001,  
                                                                                                  0.000001,  
                                                                                                  0.0,  
                                                                                                  0.0,  
                                                                                                  0.0,
                                                                                                  0.0,
                                                                                                  0.0,  
                                                                                                  0.0);
	leptonFSRCovMatrix.push_back(photonCovarianceTMatrixD);
	leptonHepFSRCovMatrix.push_back(photonCovarianceHepTMatrixD);

	leptonFSRCovMatrixID.push_back(photonCovarianceTMatrixD);
	leptonHepFSRCovMatrixID.push_back(photonCovarianceHepTMatrixD);

	leptonFSRCovMatrixMS.push_back(photonCovarianceTMatrixD);
	leptonHepFSRCovMatrixMS.push_back(photonCovarianceHepTMatrixD);

	if(isZ1)
	{
		getZ1()->setHasFSR(true);
		getZ1()->setFSR4Momentum(momFSR);
		getZ1()->setFSRError(photonCovarianceTMatrixD);
		getZ1()->setHepFSRError(photonCovarianceHepTMatrixD);
	}
	else if(isZ2)
	{
		getZ2()->setHasFSR(true);
		getZ2()->setFSR4Momentum(momFSR);
		getZ2()->setFSRError(photonCovarianceTMatrixD);
		getZ2()->setHepFSRError(photonCovarianceHepTMatrixD);
	}

}


