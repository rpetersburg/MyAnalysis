#include <stdlib.h>
#include <string>
#include "MyAnalysis/DiLepton.h"
#include <iostream>


using namespace std;
////////////////////////////////////////////////////////////////////////////////////////
//							Constructor and Destructor
////////////////////////////////////////////////////////////////////////////////////////

DiLepton::DiLepton()
{
	Reset();
}

DiLepton::DiLepton(ChargedLepton *lepOne, ChargedLepton *lepTwo)
{
	Reset();
	Set(lepOne, lepTwo);
}
DiLepton::~DiLepton()
{
	delete m_momentum_main;
	leptonInfo.clear();
	leptonLorentz.clear();		
	leptonCovMatrix.clear();
	leptonHepCovMatrix.clear();
	leptonZMassFSRLorentz.clear();
	leptonZMassCovMatrix.clear();
}
////////////////////////////////////////////////////////////////////////////////////////
//									Helpers				
////////////////////////////////////////////////////////////////////////////////////////
<<<<<<< HEAD
// Resets all the variables
=======
// Resets all the varibles
>>>>>>> bc7b9ddaf72f0a41dfe1bb5d9068cc4b03444c0d
void DiLepton::Reset()
{
	flavor = -1; 
	neutral = -1; 

	type = -1;
	
	mass = -999;
	massErr = -999;

	m_momentum_main = new TLorentzVector();

	lepPlus = 0;
	lepNeg = 0;

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

	hasFSR = false;
}

// Sets all the variable and ensures an order to the system
void DiLepton::Set(ChargedLepton *lepOne, ChargedLepton *lepTwo)
{
	flavor = lepOne->getFlavor();
	neutral = (lepOne->getCharge() * lepTwo->getCharge() < 0);

	*m_momentum_main = *lepOne->get4Momentum() + *lepTwo->get4Momentum();
	m_momentum = lepOne->get4MomentumNoP() + lepTwo->get4MomentumNoP();

	if(lepOne->getCharge() >= 0)
	{
		lepPlus = lepOne;
		lepNeg = lepTwo;
	}
	else
	{
		lepPlus = lepTwo;
		lepNeg = lepOne;
	}
	leptonInfo.push_back(getLepPlus());	
	leptonInfo.push_back(getLepNeg());	

	leptonLorentz.push_back(getLepPlus()->get4MomentumNoP());
	leptonLorentz.push_back(getLepNeg()->get4MomentumNoP());
}
// Compares if any of the leptons are same in the two class
Bool_t DiLepton::IsOverlap (DiLepton *toCompare)
{
	Bool_t testOverlap = false;
	if((getLepPlus() == toCompare->getLepPlus()) || 
			(getLepNeg() == toCompare->getLepNeg()) || 
			(getLepNeg() == toCompare->getLepPlus()) || 
			(getLepPlus() == toCompare->getLepNeg()))
	{testOverlap = true;}
	return testOverlap;
}


void DiLepton::SetElRescale(AtlasRoot::egammaEnergyCorrectionTool *telRescale)
{
	elRescale = telRescale;
	getLepPlus()->SetElRescale(telRescale);
	getLepNeg()->SetElRescale(telRescale);
}

void DiLepton::FillCovMatrix(Int_t runNumber_sf)
{
	getLepPlus()->FillCovMatrix(runNumber_sf);
	getLepNeg()->FillCovMatrix(runNumber_sf);
	
	// Order is important
	leptonCovMatrix.push_back(getLepPlus()->getCovMatrix());
	leptonCovMatrix.push_back(getLepNeg()->getCovMatrix());

	leptonHepCovMatrix.push_back(getLepPlus()->getHepCovMatrix());
	leptonHepCovMatrix.push_back(getLepNeg()->getHepCovMatrix());

}
