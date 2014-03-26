#include <stdlib.h>
#include <string>
#include "MyAnalysis/Correction.h"
#include <iostream>


using namespace std;
////////////////////////////////////////////////////////////////////////////////////////
//							Constructor and Destructor
////////////////////////////////////////////////////////////////////////////////////////

Correction::Correction(Int_t tdataYear,Bool_t tIsMC, Bool_t tSmearMC, Bool_t tScaleData, 
		Bool_t tScaleCrack, Bool_t tJetCalibration, Bool_t tScaleMuon, Bool_t tscaleEfficiency, Bool_t tuseMoriond,
		Int_t tcurrCollection, Int_t tcurrDataCollection, Bool_t tepCombination, Root::TPileupReweighting* tpileupTool)
{
	dataYear = tdataYear;
	isMC = tIsMC;
	currCollection = tcurrCollection;
	currDataCollection = tcurrDataCollection;
	smearMC = tSmearMC;
	scaleData = tScaleData;
	scaleCrack = tScaleCrack;
	jetCalibration = tJetCalibration;
	scaleMuon = tScaleMuon;
	scaleEfficiency = tscaleEfficiency;
	useMoriond = tuseMoriond;
	isDebugCall = false;
	pileupTool = tpileupTool;
	epCombination = tepCombination;
	// Init it to zero for consistency
	muSmear = 0;
	mResolSF = 0;
	elRescale = 0;
	calibrationJES = 0;
	egSFClassID = 0;
	egSFClassReco = 0;
	StacoSCF = 0;
	StacoSASCF = 0;
	CaloMuSCF = 0;

	// To print
	printInfoMuSmear = 0;
	printInfoElSmear = 0;	
	printInfoPhSmear = 0;	
	printInfoD0Z0Smear = 0;
	printInfoJetCal = 0;
	printInfoElEff = 0;
	printInfoMuEff = 0;

}
Correction::~Correction()
{
	if(muSmear != 0) delete muSmear;
	if(elRescale != 0) delete elRescale;
	if(calibrationJES != 0) delete calibrationJES;
	if(mResolSF != 0) delete mResolSF;
	if(egSFClassID != 0) delete egSFClassID;
	if(egSFClassReco != 0) delete egSFClassReco;	
	if(StacoSCF != 0) delete StacoSCF;	
	if(StacoSASCF != 0) delete StacoSASCF;	
	if(CaloMuSCF != 0) delete CaloMuSCF;	
	
}
////////////////////////////////////////////////////////////////////////////////////////
//							Muon smearing
////////////////////////////////////////////////////////////////////////////////////////
void Correction::InitMuonSmear()
{
	// For Muon Smearing
	if(dataYear == 2011)
	{
        muSmear = new MuonSmear::SmearingClass("Data11","staco","q_pT","Rel17","../../MuonMomentumCorrections/share/");
	}
	else if(dataYear == 2012)
	{
		//muSmear = new MuonSmear::SmearingClass("Data12","staco","q_pT","Rel17.2Repro","../../MuonMomentumCorrections/share/");
		muSmear = new MuonSmear::SmearingClass("Data12","staco","q_pT","Rel17.2Sum13","../../MuonMomentumCorrections/share/");
	}
	else
	{
		cout<<"Error Correction::InitMuonSmear: dataYear not reconginized"<<endl;
		return;
	}
	muSmear->UseScale(1);
	muSmear->UseImprovedCombine();
	// For Muon error
	if(!isMC)
	{
    	if(dataYear == 2012){
      	mResolSF=new Analysis::MuonResolutionAndMomentumScaleFactors("../../MuonMomentumCorrections/share/final_scale_factors_data2012.txt");
    	}
   	 	if(dataYear == 2011){
      	mResolSF=new Analysis::MuonResolutionAndMomentumScaleFactors("../../MuonMomentumCorrections/share/final_scale_factors_data2011.txt");
    	}
  	}
	else
	{
    	if(dataYear == 2012){
      	mResolSF=new Analysis::MuonResolutionAndMomentumScaleFactors("../../MuonMomentumCorrections/share/final_scale_factors_MC12_smearing.txt");
    	}
    	if(dataYear == 2011){
      	mResolSF=new Analysis::MuonResolutionAndMomentumScaleFactors("../../MuonMomentumCorrections/share/final_scale_factors_MC11_smeared.txt");
    	}
 	}
	// For Muon Efficiency... Code from Fabien
	Analysis::AnalysisMuonConfigurableScaleFactors::Configuration config;
	Analysis::AnalysisMuonConfigurableScaleFactors::Configuration config_sa;
	if(dataYear == 2011)
	{
	  config = Analysis::AnalysisMuonConfigurableScaleFactors::AverageOverPeriods;
	  config_sa = Analysis::AnalysisMuonConfigurableScaleFactors::Default;
	}
	else if (dataYear == 2012)
	{
	  //config = Analysis::AnalysisMuonConfigurableScaleFactors::AverageOverRuns;
	  //config_sa = Analysis::AnalysisMuonConfigurableScaleFactors::AverageOverRuns;
		config=Analysis::AnalysisMuonConfigurableScaleFactors::AverageOverPeriods;
		config_sa=Analysis::AnalysisMuonConfigurableScaleFactors::AverageOverPeriods;

	}
  	string directory("../../MuonEfficiencyCorrections/share/");      
	string unit("MeV");

	if(dataYear == 2011)
	{
		string muon_type_staco   ("STACO_CB_plus_ST_2011_SF.txt.gz");       // for STACO combined muon
		string muon_type_staco_sa("STACOHighEta.txt.gz");                   // for STACO combined muon
		string muon_type_calo_mu ("CaloTag_2011_SF.txt.gz");                // CaloTag (Calo muons)
		// Getting the period Luminosity to give to the muon scale factors
		Double_t periodB =	pileupTool->GetIntegratedLumi(177986,178109);
		Double_t periodD =	pileupTool->GetIntegratedLumi(179710,180481);
		Double_t periodE =	pileupTool->GetIntegratedLumi(180614,180776);
		Double_t periodF =	pileupTool->GetIntegratedLumi(182013,182519);
		Double_t periodG =	pileupTool->GetIntegratedLumi(182726,183462);
		Double_t periodH =	pileupTool->GetIntegratedLumi(183544,184169);
		Double_t periodI =	pileupTool->GetIntegratedLumi(185353,186493);
		Double_t periodJ =	pileupTool->GetIntegratedLumi(186516,186755);
		Double_t periodK =	pileupTool->GetIntegratedLumi(186873,187815);
		Double_t periodL =	pileupTool->GetIntegratedLumi(188902,190343);
		Double_t periodM =	pileupTool->GetIntegratedLumi(190503,191933);
		Double_t periodAll= pileupTool->GetIntegratedLumi();
        
		StacoSCF = new Analysis::AnalysisMuonConfigurableScaleFactors  (directory,muon_type_staco,unit,config);
        StacoSCF->addPeriod("B", periodB);
        StacoSCF->addPeriod("D", periodD);
        StacoSCF->addPeriod("E", periodE);
        StacoSCF->addPeriod("F", periodF);
        StacoSCF->addPeriod("G", periodG);
        StacoSCF->addPeriod("H", periodH);
        StacoSCF->addPeriod("I", periodI);
        StacoSCF->addPeriod("J", periodJ);
        StacoSCF->addPeriod("K", periodK);
        StacoSCF->addPeriod("L", periodL);
        StacoSCF->addPeriod("M", periodM);
		StacoSCF->Initialise();	

        StacoSASCF = new Analysis::AnalysisMuonConfigurableScaleFactors  (directory,muon_type_staco_sa,unit,config_sa); 
		StacoSASCF->addPeriod("B", periodB);
        StacoSASCF->addPeriod("D", periodD);
        StacoSASCF->addPeriod("E", periodE);
        StacoSASCF->addPeriod("F", periodF);
        StacoSASCF->addPeriod("G", periodG);
        StacoSASCF->addPeriod("H", periodH);
        StacoSASCF->addPeriod("I", periodI);
        StacoSASCF->addPeriod("J", periodJ);
        StacoSASCF->addPeriod("K", periodK);
        StacoSASCF->addPeriod("L", periodL);
        StacoSASCF->addPeriod("M", periodM);
		StacoSASCF->Initialise();

        CaloMuSCF = new Analysis::AnalysisMuonConfigurableScaleFactors (directory,muon_type_calo_mu,unit,config); 
		CaloMuSCF->addPeriod("B", periodB);
        CaloMuSCF->addPeriod("D", periodD);
        CaloMuSCF->addPeriod("E", periodE);
        CaloMuSCF->addPeriod("F", periodF);
        CaloMuSCF->addPeriod("G", periodG);
        CaloMuSCF->addPeriod("H", periodH);
        CaloMuSCF->addPeriod("I", periodI);
        CaloMuSCF->addPeriod("J", periodJ);
        CaloMuSCF->addPeriod("K", periodK);
        CaloMuSCF->addPeriod("L", periodL);
        CaloMuSCF->addPeriod("M", periodM);
		CaloMuSCF->Initialise();
	}
	else if(dataYear == 2012)
	{
		string muon_type_staco    ("STACO_CB_plus_ST_2012_SF.txt.gz");       // for STACO combined muon
		string muon_type_staco_sa ("STACO_CB_plus_ST_2012_SFms.txt.gz");     // for STACO combined muon
		string muon_type_calo_mu  ("CaloTag_2012_SF.txt.gz");                // CaloTag (Calo muons)

		Double_t period2012A	= pileupTool->GetIntegratedLumi(200804,201556);
      	Double_t period2012B	= pileupTool->GetIntegratedLumi(202660,205113);
      	Double_t period2012C	= pileupTool->GetIntegratedLumi(206248,207397);
      	Double_t period2012D	= pileupTool->GetIntegratedLumi(207447,209025);
      	Double_t period2012E	= pileupTool->GetIntegratedLumi(209074,210308);
      	Double_t period2012F	= 0.0;
      	Double_t period2012G	= pileupTool->GetIntegratedLumi(211522,212272);
      	Double_t period2012H	= pileupTool->GetIntegratedLumi(212619,213359);
      	Double_t period2012I	= pileupTool->GetIntegratedLumi(213431,213819);
      	Double_t period2012J	= pileupTool->GetIntegratedLumi(213900,215091); //
      	Double_t period2012K	= 0.0;
      	Double_t period2012L	= pileupTool->GetIntegratedLumi(215414,215643);
      	Double_t period2012M	= pileupTool->GetIntegratedLumi(216399,216432); // 2012 period A-M
      	Double_t period2012All	= pileupTool->GetIntegratedLumi();

        StacoSCF = new Analysis::AnalysisMuonConfigurableScaleFactors  (directory,muon_type_staco,unit,config);
		StacoSCF->addPeriod("A", period2012A);
        StacoSCF->addPeriod("B", period2012B);
        StacoSCF->addPeriod("C", period2012C);
        StacoSCF->addPeriod("D", period2012D);
        StacoSCF->addPeriod("E", period2012E);
        StacoSCF->addPeriod("G", period2012G);
        StacoSCF->addPeriod("H", period2012H);
        StacoSCF->addPeriod("I", period2012I);
        StacoSCF->addPeriod("J", period2012J);
        StacoSCF->addPeriod("L", period2012L);
        StacoSCF->Initialise();		
        
		StacoSASCF = new Analysis::AnalysisMuonConfigurableScaleFactors  (directory,muon_type_staco_sa,unit,config_sa);
		StacoSASCF->addPeriod("A", period2012A);
        StacoSASCF->addPeriod("B", period2012B);
        StacoSASCF->addPeriod("C", period2012C);
        StacoSASCF->addPeriod("D", period2012D);
        StacoSASCF->addPeriod("E", period2012E);
        StacoSASCF->addPeriod("G", period2012G);
        StacoSASCF->addPeriod("H", period2012H);
        StacoSASCF->addPeriod("I", period2012I);
        StacoSASCF->addPeriod("J", period2012J);
        StacoSASCF->addPeriod("L", period2012L);
		StacoSASCF->Initialise();		
        
		CaloMuSCF = new Analysis::AnalysisMuonConfigurableScaleFactors (directory,muon_type_calo_mu,unit,config); 
		CaloMuSCF->addPeriod("A", period2012A);
        CaloMuSCF->addPeriod("B", period2012B);
        CaloMuSCF->addPeriod("C", period2012C);
        CaloMuSCF->addPeriod("D", period2012D);
        CaloMuSCF->addPeriod("E", period2012E);
        CaloMuSCF->addPeriod("G", period2012G);
        CaloMuSCF->addPeriod("H", period2012H);
        CaloMuSCF->addPeriod("I", period2012I);
        CaloMuSCF->addPeriod("J", period2012J);
        CaloMuSCF->addPeriod("L", period2012L);
		CaloMuSCF->Initialise();
	}
}
void Correction::InitElectronSmear(Int_t electronCollection)
{
	// Smearing function
	//elRescale = new egRescaler::EnergyRescalerUpgrade();
	elRescale= new AtlasRoot::egammaEnergyCorrectionTool();

	if(dataYear == 2011){
	  // New Rescaler to do a lot of things related to electron/scale crack and other things
      elRescale->setFileName("../../ElectronPhotonFourMomentumCorrection/data/egammaEnergyCorrectionData.root");

	  if (currCollection == MCCollection::MC11c || currDataCollection == dataCalibType::y2011c) 
	  {
		  elRescale->setESModel(egEnergyCorr::es2011c);
	  }
	  else if (currCollection == MCCollection::MC11d || currDataCollection == dataCalibType::y2011d) 
	  {
		  elRescale->setESModel(egEnergyCorr::es2011d);
	  }
	  elRescale->initialize();
    }
	else if (dataYear == 2012 ){
	  // New Rescaler to do a lot of things related to electron/scale crack and other things
      elRescale->setFileName("../../ElectronPhotonFourMomentumCorrection/data/egammaEnergyCorrectionData.root");

      if (currCollection == MCCollection::MC12a) 
	  {
		  elRescale->setESModel(egEnergyCorr::es2012a);
		  //elRescale->useIntermoduleCorrection(false);
		  //elRescale->usePhiUniformCorrection(false);
	  }
	  else if(currCollection == MCCollection::MC12b)
	  {
		  elRescale->setESModel(egEnergyCorr::es2012a);
	  }
	  else if(currCollection == MCCollection::MC12c)
	  {
		  elRescale->setESModel(egEnergyCorr::es2012c);
	  }
	  else if (currDataCollection == dataCalibType::y2012ab)
	  {
		  elRescale->setESModel(egEnergyCorr::es2012a);
		  cout<<"Data collection: Setting flags for elRescale"<<endl;
		  elRescale->useIntermoduleCorrection(false);
		  elRescale->usePhiUniformCorrection(false);
		  elRescale->useGainCorrection(false);
	  }
	  else if (currDataCollection == dataCalibType::y2012c)
	  {
		  elRescale->setESModel(egEnergyCorr::es2012c);
		  cout<<"Data collection: Setting flags for elRescale"<<endl;		  
		  elRescale->useIntermoduleCorrection(false);
		  elRescale->usePhiUniformCorrection(false);
		  elRescale->useGainCorrection(false);
	  }

	  else cout<<"Error Correction::InitElectronSmear: currCollection not reconginized"<<endl;

	  elRescale->initialize();
    }
	else {cout<<"Error Correction::InitElectronSmear: dataYear not reconginized"<<endl;}

	// ID and reco effciency
	egSFClassID = new Root::TElectronEfficiencyCorrectionTool();
    egSFClassReco = new Root::TElectronEfficiencyCorrectionTool();

	if(dataYear == 2011)
	{
		// Reco
		if(useMoriond)
		{	
			cout<<"Electron SF config file: efficiencySF.offline.RecoTrk.2011.7TeV.rel17p0.v01.root"<<endl;
        	egSFClassReco->addFileName("../../ElectronEfficiencyCorrection/data/efficiencySF.offline.RecoTrk.2011.7TeV.rel17p0.v01.root");
		}
		else
		{
			cout<<"Electron SF config file: efficiencySF.offline.RecoTrk.2011.7TeV.rel17p0.v03.root"<<endl;
        	egSFClassReco->addFileName("../../ElectronEfficiencyCorrection/data/efficiencySF.offline.RecoTrk.2011.7TeV.rel17p0.v03.root");
			//cout<<"Electron SF config file: efficiencySF.offline.RecoTrk.2011.7TeV.rel17p0.v01.root"<<endl;
        	//egSFClassReco->addFileName("../../ElectronEfficiencyCorrection/data/efficiencySF.offline.RecoTrk.2011.7TeV.rel17p0.v01.root");

		}
		// ID
		if(electronCollection == electronCollection::LoosePlusPlus)
		{
			if(useMoriond)
			{	
				cout<<"Electron SF config file: efficiencySF.offline.Loose.2011.7TeV.rel17p0.v01.root"<<endl;
    	    	egSFClassID->addFileName("../../ElectronEfficiencyCorrection/data/efficiencySF.offline.Loose.2011.7TeV.rel17p0.v01.root");
			}
			else
			{
				cout<<"Electron SF config file: efficiencySF.offline.Loose.2011.7TeV.rel17p0.v02.root"<<endl;
    	    	egSFClassID->addFileName("../../ElectronEfficiencyCorrection/data/efficiencySF.offline.Loose.2011.7TeV.rel17p0.v02.root");
				//cout<<"Electron SF config file: efficiencySF.offline.Loose.2011.7TeV.rel17p0.v01.root"<<endl;
    	    	//egSFClassID->addFileName("../../ElectronEfficiencyCorrection/data/efficiencySF.offline.Loose.2011.7TeV.rel17p0.v01.root");

			}
		}
		else
		{
			cout<<"Electron SF config file: Not Defined"<<endl;
		}
	}
	else if(dataYear == 2012)
	{
		if(electronCollection == electronCollection::MultiLepton)
		{
			cout<<"Electron SF config file: efficiencySF.offline.Multilepton.2012.8TeV.rel17p2.v02.root"<<endl;
        	egSFClassID->addFileName("../../ElectronEfficiencyCorrection/data/efficiencySF.offline.Multilepton.2012.8TeV.rel17p2.v02.root");
			cout<<"Electron SF config file : efficiencySF.offline.RecoTrk.2012.8TeV.rel17p2.v02.root"<<endl;
        	egSFClassReco->addFileName("../../ElectronEfficiencyCorrection/data/efficiencySF.offline.RecoTrk.2012.8TeV.rel17p2.v02.root");
		}
		else if(electronCollection == electronCollection::Likelihood)
		{
			if(currCollection == MCCollection::MC12a || currCollection == MCCollection::MC12b || currDataCollection == dataCalibType::y2012ab)
			{
				cout<<"(Warning Electron) SF config file: Need to fix this"<<endl;
				cout<<"Electron SF config file: ElectronEfficiencyCorrection/data/efficiencySF.offline.LooseLLH.2012.8TeV.rel17p2.v07.root"<<endl;
        		egSFClassID->addFileName("../../ElectronEfficiencyCorrection/data/efficiencySF.offline.LooseLLH.2012.8TeV.rel17p2.v07.root");
				cout<<"Electron SF config file : ElectronEfficiencyCorrection/data/efficiencySF.offline.RecoTrk.2012.8TeV.rel17p2.GEO20.v08.root"<<endl;
        		egSFClassReco->addFileName("../../ElectronEfficiencyCorrection/data/efficiencySF.offline.RecoTrk.2012.8TeV.rel17p2.GEO20.v08.root");
			}
			else if(currCollection == MCCollection::MC12c || currDataCollection == dataCalibType::y2012c)
			{
				cout<<"Electron SF config file: ElectronEfficiencyCorrection/data/efficiencySF.offline.LooseLLH.2012.8TeV.rel17p2.GEO21.v01.root"<<endl;
        		egSFClassID->addFileName("../../ElectronEfficiencyCorrection/data/efficiencySF.offline.LooseLLH.2012.8TeV.rel17p2.GEO21.v01.root");
				cout<<"Electron SF config file : ElectronEfficiencyCorrection/data/efficiencySF.offline.RecoTrk.2012.8TeV.rel17p2.GEO21.v01.root"<<endl;
        		egSFClassReco->addFileName("../../ElectronEfficiencyCorrection/data/efficiencySF.offline.RecoTrk.2012.8TeV.rel17p2.GEO21.v01.root");
			}

		}
		else
		{
			cout<<"Electron SF config file: Not Defined"<<endl;
		}
	}

	if (!egSFClassID->initialize()) cout <<"Egamma ID SF initialization failed ..."<<endl;
    if (!egSFClassReco->initialize()) cout <<" Egamma Reco SF initialization failed ..."<<endl;

	//EP combination
    if(dataYear == 2011)
	{
      if		(currCollection == MCCollection::MC11c) 		EPcombination = new egammaFourMomentumError("2011,likelihood,MC11c");
	  else if 	(currCollection == MCCollection::MC11d) 		EPcombination = new egammaFourMomentumError("2011,likelihood,MC11d");
	  
	  else if	(currDataCollection == dataCalibType::y2011c) 	EPcombination = new egammaFourMomentumError("2011,likelihood,MC11c");
	  else if 	(currDataCollection == dataCalibType::y2011d) 	EPcombination = new egammaFourMomentumError("2011,likelihood,MC11d");

    }
	else if(dataYear == 2012)
	{
      if		(currCollection == MCCollection::MC12a) 		EPcombination = new egammaFourMomentumError("2012,likelihood,MC12a");
	  else if	(currCollection == MCCollection::MC12b) 		EPcombination = new egammaFourMomentumError("2012,likelihood,MC12a");
      else if	(currCollection == MCCollection::MC12c) 		EPcombination = new egammaFourMomentumError("2012,likelihood,MC12c");

	  else if	(currDataCollection == dataCalibType::y2012ab) 	EPcombination = new egammaFourMomentumError("2012,likelihood,MC12a");
      else if	(currDataCollection == dataCalibType::y2012c) 	EPcombination = new egammaFourMomentumError("2012,likelihood,MC12c");

    }
    
}

vector<Double_t> Correction::SmearElectron(D3PDReader::ElectronD3PDObject * el, Int_t EventNumber, Int_t RunNumber, Int_t sysVar)
{
	// just to Print a message
	if(printInfoElSmear == 0)
	{
		cout<<"--------------------------------------"<<endl;
		cout<<"El Smearing Applied"<<endl;
		cout<<"--------------------------------------"<<endl;		
		printInfoElSmear ++;
	}

	// Vector to store the smear values
	vector<Double_t> elSmearVal;
	// Clearing the Eff vector
	electronEff.clear();
	
	// Acutal Loop
	for(Int_t i = 0; i < el->n(); i++)
	{
		Double_t eta = (*el)[i].cl_eta();
		//Double_t etas2 = (*el)[i].etas2();
		Double_t etaTrk = (*el)[i].tracketa();
		//Double_t phi = (*el)[i].cl_phi();
		Double_t E = (*el)[i].cl_E();
		Double_t clPt = (*el)[i].cl_pt();
		Double_t Et = E/cosh(etaTrk);

		Double_t corrEl = 1;
		Double_t scaleCrackE = 1;

		// Saving the unsmeared E
		(*el)[i].cl_E_unsmeared = E ;

		//if(isDebugCall && sysVar == doSys::Nom)
		//{
		//	cout<<"Electron i: "<<i<<endl;
		//	cout<<"Orignal E: "<<E<<" Before Et: "<<Et<<" Before corrEl: "<<corrEl<<endl;
		//}
		// Sys
		TLorentzVector track_sys;	
		// New Rescaler tool
		if((isMC && smearMC) || (!isMC && scaleData))
		{
			PATCore::ParticleDataType::DataType dataType;
    		PATCore::ParticleType::Type ptype;
			egEnergyCorr::Scale::Variation scaleVar;
    		egEnergyCorr::Resolution::Variation res;
			egEnergyCorr::Resolution::resolutionType resType;
    		
			ptype = PATCore::ParticleType::Electron;
    		res = egEnergyCorr::Resolution::Nominal;
			scaleVar = egEnergyCorr::Scale::Nominal;
			resType = egEnergyCorr::Resolution::SigmaEff90;
		   	if(!isMC)
			{
		      	dataType=PATCore::ParticleDataType::Data; 
		    }
			else
			{
		      	dataType=PATCore::ParticleDataType::Full;
		    }
		    

		
				 if(sysVar == doSys::ZeeStatUp)			scaleVar = egEnergyCorr::Scale::ZeeStatUp;
			else if(sysVar == doSys::ZeeStatDown)		scaleVar = egEnergyCorr::Scale::ZeeStatDown;
			else if(sysVar == doSys::ZeeSystUp)			scaleVar = egEnergyCorr::Scale::ZeeSystUp;
			else if(sysVar == doSys::ZeeSystDown)		scaleVar = egEnergyCorr::Scale::ZeeSystDown;
			else if(sysVar == doSys::ZeeAllUp)			scaleVar = egEnergyCorr::Scale::ZeeAllUp;
			else if(sysVar == doSys::ZeeAllDown)		scaleVar = egEnergyCorr::Scale::ZeeAllDown;
			else if(sysVar == doSys::PSUp)				scaleVar = egEnergyCorr::Scale::PSUp;
			else if(sysVar == doSys::PSDown)			scaleVar = egEnergyCorr::Scale::PSDown;
			else if(sysVar == doSys::S12Up)				scaleVar = egEnergyCorr::Scale::S12Up;
			else if(sysVar == doSys::S12Down)			scaleVar = egEnergyCorr::Scale::S12Down;
			else if(sysVar == doSys::MatIDUp)			scaleVar = egEnergyCorr::Scale::MatIDUp;
			else if(sysVar == doSys::MatIDDown)			scaleVar = egEnergyCorr::Scale::MatIDDown;
			else if(sysVar == doSys::MatCryoUp)			scaleVar = egEnergyCorr::Scale::MatCryoUp;
			else if(sysVar == doSys::MatCryoDown)		scaleVar = egEnergyCorr::Scale::MatCryoDown;
			else if(sysVar == doSys::MatCaloUp)			scaleVar = egEnergyCorr::Scale::MatCaloUp;
			else if(sysVar == doSys::MatCaloDown)		scaleVar = egEnergyCorr::Scale::MatCaloDown;
			else if(sysVar == doSys::LArCalibUp)		scaleVar = egEnergyCorr::Scale::LArCalibUp;
			else if(sysVar == doSys::LArCalibDown)		scaleVar = egEnergyCorr::Scale::LArCalibDown;
			else if(sysVar == doSys::LArUnconvCalibUp)	scaleVar = egEnergyCorr::Scale::LArUnconvCalibUp;
			else if(sysVar == doSys::LArUnconvCalibDown)scaleVar = egEnergyCorr::Scale::LArUnconvCalibDown;
			else if(sysVar == doSys::LArElecUnconvUp)	scaleVar = egEnergyCorr::Scale::LArElecUnconvUp;
			else if(sysVar == doSys::LArElecUnconvDown)	scaleVar = egEnergyCorr::Scale::LArElecUnconvDown;
			else if(sysVar == doSys::LArElecCalibUp)	scaleVar = egEnergyCorr::Scale::LArElecCalibUp;
			else if(sysVar == doSys::LArElecCalibDown)	scaleVar = egEnergyCorr::Scale::LArElecCalibDown;
			else if(sysVar == doSys::GainUp)			scaleVar = egEnergyCorr::Scale::GainUp;
			else if(sysVar == doSys::GainDown)			scaleVar = egEnergyCorr::Scale::GainDown;
			else if(sysVar == doSys::G4Up)				scaleVar = egEnergyCorr::Scale::G4Up;
			else if(sysVar == doSys::G4Down)			scaleVar = egEnergyCorr::Scale::G4Down;
			else if(sysVar == doSys::MomentumUp)		scaleVar = egEnergyCorr::Scale::MomentumUp;
			else if(sysVar == doSys::MomentumDown)		scaleVar = egEnergyCorr::Scale::MomentumDown;
			else if(sysVar == doSys::ZSmearingUp)		res = egEnergyCorr::Resolution::ZSmearingUp;
			else if(sysVar == doSys::ZSmearingDown)		res = egEnergyCorr::Resolution::ZSmearingDown;
			else if(sysVar == doSys::SamplingTermUp)	res = egEnergyCorr::Resolution::SamplingTermUp;
			else if(sysVar == doSys::SamplingTermDown)	res = egEnergyCorr::Resolution::SamplingTermDown;
			else if(sysVar == doSys::MaterialIDUp)		res = egEnergyCorr::Resolution::MaterialIDUp;
			else if(sysVar == doSys::MaterialIDDown)	res = egEnergyCorr::Resolution::MaterialIDDown;
			else if(sysVar == doSys::MaterialCaloUp)	res = egEnergyCorr::Resolution::MaterialCaloUp;
			else if(sysVar == doSys::MaterialCaloDown)	res = egEnergyCorr::Resolution::MaterialCaloDown;
			else if(sysVar == doSys::MaterialGapUp)		res = egEnergyCorr::Resolution::MaterialGapUp;
			else if(sysVar == doSys::MaterialGapDown)	res = egEnergyCorr::Resolution::MaterialGapDown;
			else if(sysVar == doSys::MaterialCryoUp)	res = egEnergyCorr::Resolution::MaterialCryoUp;
			else if(sysVar == doSys::MaterialCryoDown)	res = egEnergyCorr::Resolution::MaterialCryoDown;
			else if(sysVar == doSys::PileUpUp)			res = egEnergyCorr::Resolution::PileUpUp;
			else if(sysVar == doSys::PileUpDown)		res = egEnergyCorr::Resolution::PileUpDown;

		    elRescale->setRandomSeed(EventNumber + 100*i);
 		 	//if(i == 0 )((*el)[i].cl_etaCalo)() = -0.393008;
  			//if(i == 0 )((*el)[i].cl_phiCalo)() = 0.092456;

			double newE = E;
			if ( ((*el)[i].author() == 1) || ((*el)[i].author() == 3) ){
				newE = elRescale->getCorrectedEnergy(RunNumber,
                                       				dataType,
                                       				AtlasRoot::egammaEnergyCorrectionTool::ParticleInformation(	((*el)[i].rawcl_Es0)(),
                                                                                                 			 	((*el)[i].rawcl_Es1)(),
                                                                                                  				((*el)[i].rawcl_Es2)(),
                                                                                                  				((*el)[i].rawcl_Es3)(),
                                                                                                  				eta,
                                                                                                  				((*el)[i].cl_phi)(),
                                                                                                  				etaTrk,
                                                                                                  				((*el)[i].cl_E)(),
                                                                                                  				((*el)[i].cl_etaCalo)(),
                                                                                                  				((*el)[i].cl_phiCalo)()),
                                      				scaleVar,
                                       				res,
													resType,
                                       				1.0 );

				if(sysVar == doSys::MomentumUp || sysVar == doSys::MomentumDown)
				{ 
		   			elRescale->setRandomSeed(EventNumber + 100*i);					
    				//Set track, cluster four momenta
    				track_sys.SetPtEtaPhiE(((*el)[i].trackpt()),
    			                   ((*el)[i].tracketa()),
    			                   ((*el)[i].trackphi()),
    			                   fabs(1./((*el)[i].trackqoverp())));

					double momentum = elRescale->getCorrectedMomentum(dataType, ptype,
                                                               track_sys.P(), // momentum
                                                               track_sys.Eta(), // trk_eta
                                                               scaleVar);
              		track_sys.SetPtEtaPhiM(momentum/cosh(track_sys.Eta()),
                                 	track_sys.Eta(),
                                 	track_sys.Phi(),
                                 	pdgElMass);



				}
    		}	
			//if(isDebugCall&& sysVar == doSys::Nom)
			//{
			//	cout<<"Input to rescaleTool: runnumber: "<<RunNumber<<" rawcl_Es0: "<< ((*el)[i].rawcl_Es0)();
			//	cout<<" rawcl_Es1: "<<((*el)[i].rawcl_Es1)() <<" rawcl_Es2: "<<((*el)[i].rawcl_Es2)() <<" rawcl_Es3: "<< ((*el)[i].rawcl_Es3)()<<endl;
			//	cout<<"eta: "<<eta<<" cl_phi: "<< ((*el)[i].cl_phi)()<<" etaTrk: "<< etaTrk<< " cl_E: "<< ((*el)[i].cl_E)();
			//	cout<< " cl_etaCalo: "<<((*el)[i].cl_etaCalo)() << " cl_phiCalo: "<< ((*el)[i].cl_phiCalo)()<<endl;
			//	cout<<"RunNumber: "<<RunNumber<<" dataType: "<<dataType<<" scaleVar: "<<scaleVar<<" res: "<<res<<" resType: "<<resType<<endl;
			//}
		   	Double_t newEt = newE/cosh(etaTrk);
			Double_t ratio =  newEt/Et;
			corrEl = corrEl * ratio;
			E = newE;
			Et = newEt;
			clPt *= ratio; 

		}

		// Storing the cluster Pt for FSR search
		bfEP_cl_pt.push_back(clPt);
		bfEP_cl_Et.push_back(Et);
		electronEpErr.push_back(E * elRescale->resolution( E, ((*el)[i].cl_eta()), PATCore::ParticleType::Electron, true));
		(*el)[i].bfEP_pT = Et;	

		if(epCombination)
		{

 			//if(isDebugCall&& sysVar == doSys::Nom)
			//{
			//	cout<<"BeforeEP E: "<<E<<" Before Et: "<<Et<<" Before corrEl: "<<corrEl<<endl;
			//}
			Bool_t epDone = false;	
			Double_t newE = E;
			if ( ((*el)[i].author() == 1) || ((*el)[i].author() == 3) ){
				TLorentzVector track, cluster, combined;
    			double combined_energy_error(0.);
   
    			//Set track, cluster four momenta
    			track.SetPtEtaPhiE(((*el)[i].trackpt()),
    			                   ((*el)[i].tracketa()),
    			                   ((*el)[i].trackphi()),
    			                   fabs(1./((*el)[i].trackqoverp())));

				if(sysVar == doSys::MomentumUp || sysVar == doSys::MomentumDown) track = track_sys;
   
    			cluster.SetPtEtaPhiE(E/cosh(eta),
    			                     ((*el)[i].cl_eta()),
    			                     ((*el)[i].cl_phi()),
    			                     E);
    			Double_t el_qoverp_LM = (*el)[i].trackqoverp();
    			for (Int_t j = 0; j < ((*el)[i].refittedTrack_LMqoverp()).size(); ++j){
    			  if( ((*el)[i].refittedTrack_author()).at(j) == 4){
    			    el_qoverp_LM = ((*el)[i].refittedTrack_LMqoverp()).at(j);
    			  }
    			}
				Double_t cl_error = E * elRescale->resolution( E, ((*el)[i].cl_eta()), PATCore::ParticleType::Electron, true);
				//cout<<"cl_error: "<<cl_error<<endl;
    			//Perform four momentum combination.
    			epDone = EPcombination->buildfourmom(track,
    			                            cluster,
    			                            el_qoverp_LM,
    			                            (*el)[i].trackcov_qoverp(),
    			                            (*el)[i].charge(),
											cl_error,
											combined,
    			                            combined_energy_error);
				electronEpErr[i] = combined_energy_error;
				//cout<<"------------------- "<< epDone<<"------------------"<<endl<<endl; 
				//if(isDebugCall && sysVar == doSys::Nom)
				//{
				//	cout<<"track: pt:"<<track.Pt()<<" eta: "<<track.Eta()<<" phi: "<<track.Phi()<<" E: "<<track.E()<<endl;
 				//	cout<<"Cluster: pT: "<<cluster.Pt()<<" eta: "<<cluster.Eta()<<" phi: "<<cluster.Phi()<<" E: "<<cluster.E()<<endl;
   			//		cout<<"el_qoverp_LM: "<<el_qoverp_LM<<" trackcov_qoverp: "<<(*el)[i].trackcov_qoverp()<<endl;
    		//		cout<<"E-E_combined  = "<<fabs(E-combined.E()) <<endl;
    		//		cout<<"E_combined  = "<<combined.E() <<endl;
				//}
    			newE = combined.E();

			}

			Double_t newEt = newE/cosh(etaTrk);
			Double_t ratio =  newEt/Et;
			corrEl = corrEl * ratio;
			E = newE;
			Et = newEt;
			clPt *= ratio; 

			//if(isDebugCall && sysVar == doSys::Nom) cout<<"AfterEP E: "<<E<<" Final Et: "<<Et<<" final corrEl: "<<corrEl<<endl;

		}


		// ID and Reco Efficiency
		if(isMC && scaleEfficiency)
		{
			if(printInfoElEff == 0)
			{
				cout<<"--------------------------------------"<<endl;
				cout<<"El ID and Reco Eff Calculated"<<endl;
				cout<<"--------------------------------------"<<endl;		
				printInfoElEff ++;
			}
 			Double_t elEff = 1;
	
			Double_t clEt = bfEP_cl_Et[i];
			// to avoid Error Message			
			if(fabs(eta)<= 2.47 && clEt > 7*1000)
			{

				const Root::TResult &sf_ID = egSFClassID->calculate(PATCore::ParticleDataType::Full,RunNumber,eta,clEt);
				const Root::TResult &sf_Reco = egSFClassReco->calculate(PATCore::ParticleDataType::Full,RunNumber,eta,clEt);	

				Double_t SF_eff_reco = sf_Reco.getScaleFactor();
				Double_t SF_eff_id = sf_ID.getScaleFactor();
				elEff = SF_eff_reco * SF_eff_id;
				//if(isDebugCall)
				//{
				//	cout<<"Eff eff i: "<<i<<endl;
				//	cout<<"RunNumber: "<<RunNumber<<" eta: "<<eta<<" clEt: "<<clEt<<endl;
				//	cout<<"reco SF: "<<SF_eff_reco<<" id SF: "<<SF_eff_reco<<endl;
				//}

			}

			electronEff.push_back(elEff);			
		}
	
		// Writing the data
		if((!isMC && scaleData) || (isMC && smearMC))
		{
			if(sysVar == doSys::Nom)
			{
				(*el)[i].cl_E() = E;
				(*el)[i].cl_pt() = clPt; 
			}
				(*el)[i].E_sysVar[sysVar] = E;

			//if(sysVar == doSys::Nom)
			//{
			//	cout<<"Nom: E "<<E<<" Et"<<Et<<" cl_eta: "<<eta<<" trackEta: "<<etaTrk<<endl;
			//}
			//if(sysVar == doSys::MomentumUp)
			//{
			//	cout<<"MomentumUp: E "<<E<<" Et"<<Et<<" cl_eta: "<<eta<<" trackEta: "<<etaTrk<<endl;
			//}
			//if(sysVar == doSys::MomentumDown)
			//{
			//	cout<<"MomentumDown: E "<<E<<" Et"<<Et<<" cl_eta: "<<eta<<" trackEta: "<<etaTrk<<endl;
			//}

		}
		elSmearVal.push_back(corrEl);
	}

	return elSmearVal;
}

vector<Double_t> Correction::SmearPhoton(D3PDReader::PhotonD3PDObject * ph, Int_t EventNumber, Int_t RunNumber)
{
	// just to Print a message
	if(printInfoPhSmear == 0)
	{
		cout<<"--------------------------------------"<<endl;
		cout<<"Ph Smearing Applied"<<endl;
		cout<<"--------------------------------------"<<endl;		
		printInfoPhSmear ++;
	}

	// Vector to store the smear values
	vector<Double_t> phSmearVal;
	
	// Acutal Loop
	for(Int_t i = 0; i < ph->n(); i++)
	{
		Double_t E = (*ph)[i].E();
		Double_t E_cl = (*ph)[i].cl_E();
		Double_t Et = (*ph)[i].Et();

		Double_t corrPh = 1;
		Double_t scaleCrackE = 1;
		// Saving the unsmeared E
		(*ph)[i].cl_E_unsmeared = E_cl ;
		
		//if(isDebugCall)
		//{
		//	cout<<"Electron i: "<<i<<endl;
		//	cout<<"Orignal E: "<<E<<" Before Et: "<<Et<<" Before corrEl: "<<corrEl<<endl;
		//}
		
		// New Rescaler tool
		if((isMC && smearMC) || (!isMC && scaleData))
		{
			PATCore::ParticleDataType::DataType dataType;
    		PATCore::ParticleType::Type ptype;
			egEnergyCorr::Scale::Variation scaleVar;
    		egEnergyCorr::Resolution::Variation res;
			egEnergyCorr::Resolution::resolutionType resType;

    		ptype = PATCore::ParticleType::Photon;
    		res = egEnergyCorr::Resolution::Nominal;
			scaleVar = egEnergyCorr::Scale::Nominal;
			resType = egEnergyCorr::Resolution::SigmaEff90;
			
		   	if(!isMC)
			{
		      	dataType=PATCore::ParticleDataType::Data; 
		    }
			else
			{
		      	dataType=PATCore::ParticleDataType::Full;
		    }
		    
		    //elRescaleNew->setDebug(true);
		    //std::cout<<" runnumber = "<<  runnumber <<std::endl;
		    elRescale->setRandomSeed(EventNumber + 100*i);
			
			double newE = E;
			if ( ((*ph)[i].author() == 4) || ((*ph)[i].author() == 16) ){
				newE = elRescale->getCorrectedEnergy(RunNumber,
                                       				dataType,
													AtlasRoot::egammaEnergyCorrectionTool::ParticleInformation(((*ph)[i].rawcl_Es0)(),
                                                                                                 			   ((*ph)[i].rawcl_Es1)(),
                                                                                                  			   ((*ph)[i].rawcl_Es2)(),
                                                                                                  			   ((*ph)[i].rawcl_Es3)(),
                                                                                                  			   ((*ph)[i].cl_eta)(),
                                                                                                  			   ((*ph)[i].cl_phi)(),
                                                                                                  			   ((*ph)[i].cl_E)(),
                                                                                                  			   ((*ph)[i].cl_etaCalo)(),
                                                                                                  			   ((*ph)[i].cl_phiCalo)(),
                                                                                                  			   ((*ph)[i].ptconv)(),
                                                                                                  			   ((*ph)[i].pt1conv)(),
                                                                                                  			   ((*ph)[i].pt2conv)(),
                                                                                                  			   ((*ph)[i].convtrk1nPixHits)(),
                                                                                                  			   ((*ph)[i].convtrk1nSCTHits)(),
                                                                                                  			   ((*ph)[i].convtrk2nPixHits)(),
                                                                                                  			   ((*ph)[i].convtrk2nSCTHits)(),
                                                                                                  			   ((*ph)[i].Rconv)()),
                                      				scaleVar,
                                       				res,
													resType,													
                                       				1.0 );
    		}	
                                             
			Double_t ratio =  newE/E;
			corrPh = corrPh * ratio;
			Et = Et*ratio;
			E_cl = E_cl * ratio;

		}


		
		// Writing the data
		if((!isMC && scaleData) || (isMC && smearMC))
		{
			(*ph)[i].E() = E;
			(*ph)[i].cl_E() = E_cl;
			(*ph)[i].Et() = Et;

		}
		phSmearVal.push_back(corrPh);
	}

	return phSmearVal;
}

void Correction::SmearMuon(D3PDReader::MuonD3PDObject * mu, Int_t EventNumber, Int_t type)
{
	// just to Print a message
	if(printInfoMuSmear == 0)
	{
		cout<<"--------------------------------------"<<endl;
		cout<<"Mu Smearing Applied"<<endl;
		cout<<"--------------------------------------"<<endl;		
		printInfoMuSmear ++;
	}
	// A vector to store the smearing info
	vector<Double_t> muSmearVal;
	// to store eff values
	if(type == leptonType::MuonCalo) muonCaloEff.clear();
	else if (type == leptonType::MuonStaco) muonStacoEff.clear();

	// Actual loop
	for(Int_t i = 0; i < mu->n(); i++)
	{ 
		if(type == leptonType::MuonCalo)
		{
			Double_t phi = (*mu)[i].phi();
			Double_t E = (*mu)[i].E();
			Double_t charge = (*mu)[i].charge();
			
		//	Int_t isCombinedMu = (*mu)[i].isCombinedMuon();
		//	Int_t isStandAloneMu = (*mu)[i].isStandAloneMuon();
		//	Int_t isSegmentMu = (*mu)[i].isSegmentTaggedMuon();
			Int_t isCaloIdMu = (*mu)[i].isCaloMuonId();
			
			Double_t eta = (*mu)[i].eta();
			Double_t pT = (*mu)[i].pt();

			Double_t p = 0;
			Double_t ptMs = 0;
			
			if(isCaloIdMu) ptMs = pT;
			else ptMs = (1/fabs((*mu)[i].me_qoverp())*sin((*mu)[i].me_theta()));

		//	Double_t ptid = (1/fabs((*mu)[i].id_qoverp())*sin((*mu)[i].id_theta()));
		//	Double_t pTCombinedSmeared = 0;
		//	Double_t pTMsSmeared = 0;

			// Smear MC
			Double_t smear = 1.;
			if(isMC && smearMC)
			{
				muSmear->SetSeed(EventNumber, i);

				if(isCaloIdMu && type == leptonType::MuonCalo)
				{
					muSmear->Event(pT, eta, "ID", charge, phi);

					// Get smeard Pts
					Double_t pTID_smeared = muSmear->pTID();
					
					smear = pTID_smeared/pT;
					pT = pTID_smeared;
					p = pT * sin((*mu)[i].id_theta());
					E = E * smear;

					(*mu)[i].id_pt = pTID_smeared;
					(*mu)[i].me_pt = pTID_smeared;

				}
				else {cout<<"Correction::SmearMuon: Calo Type Mismatch"<<endl;}
			}
			// Muon Efficiency
			if(isMC && scaleEfficiency)
			{
				if(printInfoMuEff == 0)
				{
					cout<<"--------------------------------------"<<endl;
					cout<<"Muon ID and Reco Eff Calculated"<<endl;
					cout<<"--------------------------------------"<<endl;		
					printInfoMuEff ++;
				}
				
				Double_t muEff = 1;
				// Calo Muon Pt Cut
				Double_t pTCutCalo = -9999;
				if(dataYear == 2011) pTCutCalo = 5 * 1000;
				else if(dataYear == 2012) pTCutCalo = 15 * 1000;

				if(pT >= pTCutCalo && fabs(eta) <= 0.4)
				{
					TLorentzVector tlv;
					tlv.SetPtEtaPhiE(pT, eta, phi, E);
					if(isCaloIdMu && dataYear == 2011) muEff = CaloMuSCF->scaleFactor(tlv);
					else if(isCaloIdMu && dataYear == 2012) muEff = CaloMuSCF->scaleFactor(charge, tlv);
					
					else cout<<"Error: Correction::SmearMuon: Eff Calo Type Mismatch"<<endl;
					if(isDebugCall)
					{
						//cout<<"--------------"<<endl;
						//cout<<"Calo Muon: "<<" pT: "<<tlv.Pt()<<" mueff: "<<muEff<<endl;
					}
				}
				muonCaloEff.push_back(muEff);
			}
			// Scaling Data
			if(!isMC && scaleData)
			{
				Double_t scaleFactor = 1; // Just in case
				pT = pT * fabs(scaleFactor);
				p = p * fabs(scaleFactor);
				E = E * fabs(scaleFactor);
			}
			
			// OverWriting the data
			if((isMC && smearMC) || (!isMC && scaleData))
			{
				(*mu)[i].pt() = pT;
				(*mu)[i].E() = E;
			}
			muSmearVal.push_back(smear);

			// Muon Uncertainity Correction
			if(scaleMuon)
			{
				// muon type, =1 for combined muons, =2 for calorimeter and segment tagged muons, =3 for stand-alone muons 
      			Int_t type=2;
				TLorentzVector muLorentz;
      			phi = (*mu)[i].phi();
				E = (*mu)[i].E();
				eta = (*mu)[i].eta();
				pT = (*mu)[i].pt();
      			muLorentz.SetPtEtaPhiE( pT, eta, phi, E);
      
      			Double_t MuonErrSF = mResolSF->getResolutionScaleFactor(muLorentz,type);
 
      			Double_t mu_cov_qoverp_exPV       =(*mu)[i].cov_qoverp_exPV()*MuonErrSF*MuonErrSF;
      			Double_t mu_cov_d0_qoverp_exPV    =(*mu)[i].cov_d0_qoverp_exPV()*MuonErrSF;
     			Double_t mu_cov_z0_qoverp_exPV    =(*mu)[i].cov_z0_qoverp_exPV()*MuonErrSF;    
     			Double_t mu_cov_phi_qoverp_exPV   =(*mu)[i].cov_phi_qoverp_exPV()*MuonErrSF;
      			Double_t mu_cov_theta_qoverp_exPV =(*mu)[i].cov_theta_qoverp_exPV()*MuonErrSF;

      			(*mu)[i].cov_qoverp_exPV() = mu_cov_qoverp_exPV;
      			(*mu)[i].cov_d0_qoverp_exPV() = mu_cov_d0_qoverp_exPV;
      			(*mu)[i].cov_z0_qoverp_exPV()   = mu_cov_z0_qoverp_exPV;
      			(*mu)[i].cov_phi_qoverp_exPV()  = mu_cov_phi_qoverp_exPV;
      			(*mu)[i].cov_theta_qoverp_exPV()= mu_cov_theta_qoverp_exPV;

				// for ID and MS
				TLorentzVector muLorentzME;
      			muLorentzME.SetPtEtaPhiM((*mu)[i].me_pt, (*mu)[i].me_eta, (*mu)[i].me_phi(), pdgMuMass);
      
      			Double_t MuonErrSFME = mResolSF->getResolutionScaleFactor(muLorentzME,type);
 
      			Double_t mu_me_cov_qoverp_exPV       =(*mu)[i].me_cov_qoverp_exPV()*MuonErrSFME*MuonErrSFME;
      			Double_t mu_me_cov_d0_qoverp_exPV    =(*mu)[i].me_cov_d0_qoverp_exPV()*MuonErrSFME;
     			Double_t mu_me_cov_z0_qoverp_exPV    =(*mu)[i].me_cov_z0_qoverp_exPV()*MuonErrSFME;    
     			Double_t mu_me_cov_phi_qoverp_exPV   =(*mu)[i].me_cov_phi_qoverp_exPV()*MuonErrSFME;
      			Double_t mu_me_cov_theta_qoverp_exPV =(*mu)[i].me_cov_theta_qoverp_exPV()*MuonErrSFME;

      			(*mu)[i].me_cov_qoverp_exPV() 		= mu_me_cov_qoverp_exPV;
      			(*mu)[i].me_cov_d0_qoverp_exPV() 	= mu_me_cov_d0_qoverp_exPV;
      			(*mu)[i].me_cov_z0_qoverp_exPV()   	= mu_me_cov_z0_qoverp_exPV;
      			(*mu)[i].me_cov_phi_qoverp_exPV()  	= mu_me_cov_phi_qoverp_exPV;
      			(*mu)[i].me_cov_theta_qoverp_exPV()	= mu_me_cov_theta_qoverp_exPV;
				
				TLorentzVector muLorentzID;
      			muLorentzID.SetPtEtaPhiM((*mu)[i].id_pt, (*mu)[i].id_eta, (*mu)[i].id_phi(), pdgMuMass);
      
      			Double_t MuonErrSFID = mResolSF->getResolutionScaleFactor(muLorentzID,type);
 
      			Double_t mu_id_cov_qoverp_exPV       =(*mu)[i].id_cov_qoverp_exPV()*MuonErrSFID*MuonErrSFID;
      			Double_t mu_id_cov_d0_qoverp_exPV    =(*mu)[i].id_cov_d0_qoverp_exPV()*MuonErrSFID;
     			Double_t mu_id_cov_z0_qoverp_exPV    =(*mu)[i].id_cov_z0_qoverp_exPV()*MuonErrSFID;    
     			Double_t mu_id_cov_phi_qoverp_exPV   =(*mu)[i].id_cov_phi_qoverp_exPV()*MuonErrSFID;
      			Double_t mu_id_cov_theta_qoverp_exPV =(*mu)[i].id_cov_theta_qoverp_exPV()*MuonErrSFID;

      			(*mu)[i].id_cov_qoverp_exPV() 		= mu_id_cov_qoverp_exPV;
      			(*mu)[i].id_cov_d0_qoverp_exPV() 	= mu_id_cov_d0_qoverp_exPV;
      			(*mu)[i].id_cov_z0_qoverp_exPV()   	= mu_id_cov_z0_qoverp_exPV;
      			(*mu)[i].id_cov_phi_qoverp_exPV()  	= mu_id_cov_phi_qoverp_exPV;
      			(*mu)[i].id_cov_theta_qoverp_exPV()	= mu_id_cov_theta_qoverp_exPV;
				//cout<<"Calo ID eta: "<<muLorentzID.Eta()<<" MS: "<<muLorentzME.Eta()<<endl;				

			}
		}// end if for type = typeMuonCalo
		else 
		{
			Double_t phi = (*mu)[i].phi();
			Double_t E = (*mu)[i].E();
			Double_t charge = (*mu)[i].charge();
			
			Int_t isCombinedMu = (*mu)[i].isCombinedMuon();
			Int_t isStandAloneMu = (*mu)[i].isStandAloneMuon();
			Int_t isSegmentMu = (*mu)[i].isSegmentTaggedMuon();
		//	Int_t isCaloIdMu = (*mu)[i].isCaloMuonId();
			
			Double_t eta = (*mu)[i].eta();
			Double_t pT = (*mu)[i].pt();

			Double_t p = 0;
			Double_t ptMs = (1/fabs(double((*mu)[i].me_qoverp()))*sin((*mu)[i].me_theta()));
			Double_t pTId;
			if(isStandAloneMu) pTId = pT;
			else pTId = (1/fabs(double((*mu)[i].id_qoverp()))*sin((*mu)[i].id_theta())) ;

			// acutal pT Smear
			Double_t smear = 1;
			if (isMC && smearMC)
			{	
				muSmear->SetSeed(EventNumber, i);
				
				// Searing for combined muons
				if(isCombinedMu && type == leptonType::MuonStaco)
				{
					Double_t pTCBSmeared = pT;
					Double_t pTMSSmeared = pT;
					Double_t pTIDSmeared = pT;
					//if(isDebugCall)
					//{
					//	cout<<"Muon "<< i<<endl;
					//	cout<<"Before Smear MS: "<<ptMs;
					//	cout<<" ID: "<<pTId;
					//	cout<<" CB: "<<pT;
					//	cout<<" Eta: "<<eta;
					//	cout<<" IDOnly: "<<(*mu)[i].id_pt;
					//	cout<<" MSOnly: "<<(*mu)[i].me_pt<<endl;

					//}
					if(fabs(eta) <= 2.7) // to supress warings
					{

						muSmear->Event(ptMs, pTId, pT, eta, charge, phi);
						pTCBSmeared = muSmear->pTCB();
						pTMSSmeared = muSmear->pTMS();
						pTIDSmeared = muSmear->pTID();

						//cout<<"Inputs to smearing function:"<<endl;
						//cout<<"id_pt: "<<(*mu)[i].id_pt<<" id_eta: "<<(*mu)[i].id_eta<<" id_phi: "<<(*mu)[i].id_phi()<<endl;
						//cout<<"me_pt: "<<(*mu)[i].me_pt<<" me_eta: "<<(*mu)[i].me_eta<<" me_phi: "<<(*mu)[i].me_phi()<<endl;
						//if(fabs((*mu)[i].me_eta) <= 2.7)
						//{
						//	muSmear->Event((*mu)[i].me_pt, (*mu)[i].me_eta, "MS", charge, (*mu)[i].me_phi());
						//	(*mu)[i].me_pt = muSmear->pTMS();
						//}

						//if(fabs((*mu)[i].id_eta) <= 2.7)
						//{
						//	muSmear->Event((*mu)[i].id_pt, (*mu)[i].id_eta, "ID", charge, (*mu)[i].id_phi());
						//	(*mu)[i].id_pt = muSmear->pTID();
						//}

						(*mu)[i].me_pt = pTMSSmeared;
						(*mu)[i].id_pt = pTIDSmeared;
						
											
						//cout<<"ouptut to smearing function:"<<endl;
						//cout<<"id_pt: "<<(*mu)[i].id_pt<<" id_eta: "<<(*mu)[i].id_eta<<" id_phi: "<<(*mu)[i].id_phi()<<endl;
						//cout<<"me_pt: "<<(*mu)[i].me_pt<<" me_eta: "<<(*mu)[i].me_eta<<" me_phi: "<<(*mu)[i].me_phi()<<endl;


					}
					//if(isDebugCall)
					//{
					//	cout<<"After Smear MS: "<<pTMSSmeared;
					//	cout<<" ID: "<<pTIDSmeared;
					//	cout<<" CB: "<<pTCBSmeared;
					//	cout<<" Eta: "<<eta;
					//	cout<<" IDOnly: "<<(*mu)[i].id_pt;
					//	cout<<" MSOnly: "<<(*mu)[i].me_pt<<endl;
					//}

					smear = pTCBSmeared/pT;
					pT = pTCBSmeared;
					E = E*smear;
				} // End for smearing for combined muon
				else if(isStandAloneMu && type == leptonType::MuonStaco) // smearing for standalone
				{
					Double_t pTMSSmeared = pT;
					if(fabs(eta) <= 2.7 )// to supress warnings
					{
						muSmear->Event(pT, eta, "MS", charge, phi);
						pTMSSmeared = muSmear->pTMS();
					}
					smear = pTMSSmeared/pT;
					pT = pTMSSmeared;
					E = E*smear;

					(*mu)[i].id_pt = pTMSSmeared;
					(*mu)[i].me_pt = pTMSSmeared;
				} // end of smearing for muon standalone
				else if(isSegmentMu && type == leptonType::MuonStaco) // smearing for segment tagged
				{ 
					Double_t pTID_smeared = pT;
					if(fabs(eta) <= 2.7 )// to supress warnings
					{
						muSmear->Event(pT, eta, "ID", charge, phi);
						pTID_smeared = muSmear->pTID();
					}
					smear = pTID_smeared/pT;
					pT = pTID_smeared;
					E = E*smear;

					
					(*mu)[i].id_pt = pTID_smeared;
					(*mu)[i].me_pt = pTID_smeared;

				} // end of smearing for segment tagged
				else {cout<<"Correction::SmearMuon: Staco Type Mismatch"<<endl;}
				//// Smearing for ID
				//if(fabs((*mu)[i].id_eta) <= 2.7 )
				//{
				//	muSmear->Event((*mu)[i].id_pt, (*mu)[i].id_eta, "ID", charge, (*mu)[i].id_phi());
				//	(*mu)[i].id_pt = muSmear->pTID();
				//}
				//// Smearing for MS
				//if(fabs((*mu)[i].me_eta) <= 2.7 )
				//{
				//	muSmear->Event((*mu)[i].me_pt, (*mu)[i].me_eta, "MS", charge, (*mu)[i].me_phi());
				//	(*mu)[i].me_pt = muSmear->pTMS();
				//}
			}// End of actual pT Smear
			
			// Muon Efficiency
			if(isMC && scaleEfficiency)
			{
								
				Double_t muEff = 1;
				//  Muon Pt Cut
				Double_t pTCut = 5*1000;
				
				if(pT >= pTCut)
				{
					TLorentzVector tlv;
					tlv.SetPtEtaPhiE(pT, eta, phi, E);

					if(isStandAloneMu)
					{
						if(fabs(eta) >= 2.5 && fabs(eta) <= 2.7)
						{
							if(dataYear == 2011) muEff = StacoSASCF->scaleFactor(tlv);
							else if(dataYear == 2012) muEff = StacoSASCF->scaleFactor(charge, tlv);
						}
					}
					else
					{
						muEff = StacoSCF->scaleFactor(charge, tlv);
					}
					//if(isDebugCall)
					//{
					//	cout<<"Staco Muon: "<<" pT: "<<tlv.Pt()<<" eta: "<<tlv.Eta()<<" phi: "<<tlv.Phi()<<" E: "<<tlv.E()<<" mueff: "<<muEff<<endl;
					//}

				}

				muonStacoEff.push_back(muEff);
			}
			// Scaling Data
			if(!isMC && scaleData)
			{
				Double_t scaleFactor = 1; // Just in case
				pT = pT * fabs(scaleFactor);
				p = p * fabs(scaleFactor);
				E = E * fabs(scaleFactor);
			}
			
			// OverWriting the data
			if((isMC && smearMC) || (!isMC && scaleData))
			{
				(*mu)[i].pt() = pT;
				(*mu)[i].E() = E;
			}
			muSmearVal.push_back(smear);			
			if (pT<1e-3  || pT>1e7)
			{
     			cout<<" WARNING : Muon Smear values are unphysical"<<endl;
    		}
			// Muon Uncertainity Correction
			if(scaleMuon)
			{
				// muon type, =1 for combined muons, =2 for calorimeter and segment tagged muons, =3 for stand-alone muons 
      			Int_t type=0;
				if (isStandAloneMu){    type = 3;}
      			else if (isSegmentMu){  type = 2;}
      			else if (isCombinedMu){ type = 1;}
				TLorentzVector muLorentz;
      			phi = (*mu)[i].phi();
				E = (*mu)[i].E();
				eta = (*mu)[i].eta();
				pT = (*mu)[i].pt();
      			muLorentz.SetPtEtaPhiE( pT, eta, phi, E);
     
				//cout<<"-----------------"<<endl;
				//cout<<std::scientific<<"pT "<<pT<<endl<<"beforeSmear"<<endl;

      			Double_t MuonErrSF = mResolSF->getResolutionScaleFactor(muLorentz,type);

				//cout<<"MuonErrSF: "<<MuonErrSF<<endl;
 				//cout<<"mu_cov_qoverp_exPV: "		<<(*mu)[i].cov_qoverp_exPV()<<endl;
 				//cout<<"mu_cov_d0_qoverp_exPV: "		<<(*mu)[i].cov_d0_qoverp_exPV()<<endl;
 				//cout<<"mu_cov_z0_qoverp_exPV: "		<<(*mu)[i].cov_z0_qoverp_exPV()<<endl;
 				//cout<<"mu_cov_phi_qoverp_exPV: "	<<(*mu)[i].cov_phi_qoverp_exPV()<<endl;
 				//cout<<"mu_cov_theta_qoverp_exPV: "	<<(*mu)[i].cov_qoverp_exPV()<<endl;

				Double_t mu_cov_qoverp_exPV       =(*mu)[i].cov_qoverp_exPV()*MuonErrSF*MuonErrSF;
      			Double_t mu_cov_d0_qoverp_exPV    =(*mu)[i].cov_d0_qoverp_exPV()*MuonErrSF;
     			Double_t mu_cov_z0_qoverp_exPV    =(*mu)[i].cov_z0_qoverp_exPV()*MuonErrSF;    
     			Double_t mu_cov_phi_qoverp_exPV   =(*mu)[i].cov_phi_qoverp_exPV()*MuonErrSF;
      			Double_t mu_cov_theta_qoverp_exPV =(*mu)[i].cov_theta_qoverp_exPV()*MuonErrSF;

      			(*mu)[i].cov_qoverp_exPV() 		= mu_cov_qoverp_exPV;
      			(*mu)[i].cov_d0_qoverp_exPV() 	= mu_cov_d0_qoverp_exPV;
      			(*mu)[i].cov_z0_qoverp_exPV()   = mu_cov_z0_qoverp_exPV;
      			(*mu)[i].cov_phi_qoverp_exPV()  = mu_cov_phi_qoverp_exPV;
      			(*mu)[i].cov_theta_qoverp_exPV()= mu_cov_theta_qoverp_exPV;

				//cout<<"afterSmear"<<endl;
 				//cout<<"mu_cov_qoverp_exPV: "		<<(*mu)[i].cov_qoverp_exPV()<<endl;
 				//cout<<"mu_cov_d0_qoverp_exPV: "		<<(*mu)[i].cov_d0_qoverp_exPV()<<endl;
 				//cout<<"mu_cov_z0_qoverp_exPV: "		<<(*mu)[i].cov_z0_qoverp_exPV()<<endl;
 				//cout<<"mu_cov_phi_qoverp_exPV: "	<<(*mu)[i].cov_phi_qoverp_exPV()<<endl;
 				//cout<<"mu_cov_theta_qoverp_exPV: "	<<(*mu)[i].cov_qoverp_exPV()<<endl;


				// for ID and MS
				TLorentzVector muLorentzME;
				muLorentzME.SetPtEtaPhiM((*mu)[i].me_pt, (*mu)[i].me_eta, (*mu)[i].me_phi(), pdgMuMass);
      
      			Double_t MuonErrSFME = mResolSF->getResolutionScaleFactor(muLorentzME,3);
      			//Double_t MuonErrSFME = mResolSF->getResolutionScaleFactor(muLorentzME,type);
 
      			Double_t mu_me_cov_qoverp_exPV       =(*mu)[i].me_cov_qoverp_exPV()*MuonErrSFME*MuonErrSFME;
      			Double_t mu_me_cov_d0_qoverp_exPV    =(*mu)[i].me_cov_d0_qoverp_exPV()*MuonErrSFME;
     			Double_t mu_me_cov_z0_qoverp_exPV    =(*mu)[i].me_cov_z0_qoverp_exPV()*MuonErrSFME;    
     			Double_t mu_me_cov_phi_qoverp_exPV   =(*mu)[i].me_cov_phi_qoverp_exPV()*MuonErrSFME;
      			Double_t mu_me_cov_theta_qoverp_exPV =(*mu)[i].me_cov_theta_qoverp_exPV()*MuonErrSFME;
				//cout<<"MuonErrSFME: "<<MuonErrSFME<<endl;
 				//cout<<"me_mu_cov_qoverp_exPV: "			<<(*mu)[i].me_cov_qoverp_exPV()<<endl;
 				//cout<<"me_mu_cov_d0_qoverp_exPV: "		<<(*mu)[i].me_cov_d0_qoverp_exPV()<<endl;
 				//cout<<"me_mu_cov_z0_qoverp_exPV: "		<<(*mu)[i].me_cov_z0_qoverp_exPV()<<endl;
 				//cout<<"me_mu_cov_phi_qoverp_exPV: "		<<(*mu)[i].me_cov_phi_qoverp_exPV()<<endl;
 				//cout<<"me_mu_cov_theta_qoverp_exPV: "	<<(*mu)[i].me_cov_qoverp_exPV()<<endl;

      			(*mu)[i].me_cov_qoverp_exPV() 		= mu_me_cov_qoverp_exPV;
      			(*mu)[i].me_cov_d0_qoverp_exPV() 	= mu_me_cov_d0_qoverp_exPV;
      			(*mu)[i].me_cov_z0_qoverp_exPV()   	= mu_me_cov_z0_qoverp_exPV;
      			(*mu)[i].me_cov_phi_qoverp_exPV()  	= mu_me_cov_phi_qoverp_exPV;
      			(*mu)[i].me_cov_theta_qoverp_exPV()	= mu_me_cov_theta_qoverp_exPV;

				//cout<<"afterSmear"<<endl;
 				//cout<<"me_mu_cov_qoverp_exPV: "			<<(*mu)[i].me_cov_qoverp_exPV()<<endl;
 				//cout<<"me_mu_cov_d0_qoverp_exPV: "		<<(*mu)[i].me_cov_d0_qoverp_exPV()<<endl;
 				//cout<<"me_mu_cov_z0_qoverp_exPV: "		<<(*mu)[i].me_cov_z0_qoverp_exPV()<<endl;
 				//cout<<"me_mu_cov_phi_qoverp_exPV: "		<<(*mu)[i].me_cov_phi_qoverp_exPV()<<endl;
 				//cout<<"me_mu_cov_theta_qoverp_exPV: "	<<(*mu)[i].me_cov_qoverp_exPV()<<endl;
				
				TLorentzVector muLorentzID;
      			muLorentzID.SetPtEtaPhiM((*mu)[i].id_pt, (*mu)[i].id_eta, (*mu)[i].id_phi(), pdgMuMass);
      
      			Double_t MuonErrSFID = mResolSF->getResolutionScaleFactor(muLorentzID,2);

				//cout<<"MuonErrSFID: "<<MuonErrSFID<<endl;
 				//cout<<"id_mu_cov_qoverp_exPV: "			<<(*mu)[i].id_cov_qoverp_exPV()<<endl;
 				//cout<<"id_mu_cov_d0_qoverp_exPV: "		<<(*mu)[i].id_cov_d0_qoverp_exPV()<<endl;
 				//cout<<"id_mu_cov_z0_qoverp_exPV: "		<<(*mu)[i].id_cov_z0_qoverp_exPV()<<endl;
 				//cout<<"id_mu_cov_phi_qoverp_exPV: "		<<(*mu)[i].id_cov_phi_qoverp_exPV()<<endl;
 				//cout<<"id_mu_cov_theta_qoverp_exPV: "	<<(*mu)[i].id_cov_qoverp_exPV()<<endl;

      			Double_t mu_id_cov_qoverp_exPV       =(*mu)[i].id_cov_qoverp_exPV()*MuonErrSFID*MuonErrSFID;
      			Double_t mu_id_cov_d0_qoverp_exPV    =(*mu)[i].id_cov_d0_qoverp_exPV()*MuonErrSFID;
     			Double_t mu_id_cov_z0_qoverp_exPV    =(*mu)[i].id_cov_z0_qoverp_exPV()*MuonErrSFID;    
     			Double_t mu_id_cov_phi_qoverp_exPV   =(*mu)[i].id_cov_phi_qoverp_exPV()*MuonErrSFID;
      			Double_t mu_id_cov_theta_qoverp_exPV =(*mu)[i].id_cov_theta_qoverp_exPV()*MuonErrSFID;
      			
				(*mu)[i].id_cov_qoverp_exPV() 		= mu_id_cov_qoverp_exPV;
      			(*mu)[i].id_cov_d0_qoverp_exPV() 	= mu_id_cov_d0_qoverp_exPV;
      			(*mu)[i].id_cov_z0_qoverp_exPV()   	= mu_id_cov_z0_qoverp_exPV;
      			(*mu)[i].id_cov_phi_qoverp_exPV()  	= mu_id_cov_phi_qoverp_exPV;
      			(*mu)[i].id_cov_theta_qoverp_exPV()	= mu_id_cov_theta_qoverp_exPV;
				//cout<<"staco ID eta: "<<muLorentzID.Eta()<<" MS: "<<muLorentzME.Eta()<<endl;	
				//cout<<"afterSmear"<<endl;
 				//cout<<"id_mu_cov_qoverp_exPV: "			<<(*mu)[i].id_cov_qoverp_exPV()<<endl;
 				//cout<<"id_mu_cov_d0_qoverp_exPV: "		<<(*mu)[i].id_cov_d0_qoverp_exPV()<<endl;
 				//cout<<"id_mu_cov_z0_qoverp_exPV: "		<<(*mu)[i].id_cov_z0_qoverp_exPV()<<endl;
 				//cout<<"id_mu_cov_phi_qoverp_exPV: "		<<(*mu)[i].id_cov_phi_qoverp_exPV()<<endl;
 				//cout<<"id_mu_cov_theta_qoverp_exPV: "	<<(*mu)[i].id_cov_qoverp_exPV()<<endl;


			}

		}// end if for type != typeMuonCalo
	}
}
////////////////////////////////////////////////////////////////////////////////////////
//							Jet cal
////////////////////////////////////////////////////////////////////////////////////////
void Correction::InitJetCal()
{
	TString jetAlgo = "AntiKt4TopoEM";
 	TString JES_config_file = "";	
	if(dataYear == 2011)
	{
    	JES_config_file="../../ApplyJetCalibration/data/CalibrationConfigs/InsituJES_2011_Preliminary.config";
    	calibrationJES = new JetAnalysisCalib::JetCalibrationTool(jetAlgo,JES_config_file, !isMC);
  	}
	else if(dataYear == 2012)
	{
		if(!isMC) 	 
			JES_config_file = "../../ApplyJetCalibration/data/CalibrationConfigs/JES_Full2012dataset_Preliminary_Jan13.config";
		else if	(currCollection == MCCollection::MC12a) 
			JES_config_file = "../../ApplyJetCalibration/data/CalibrationConfigs/JES_Full2012dataset_Preliminary_Jan13.config";
		else if	(currCollection == MCCollection::MC12b) 
			JES_config_file = "../../ApplyJetCalibration/data/CalibrationConfigs/JES_Full2012dataset_Preliminary_MC12b_Sep23.config";
		else if	(currCollection == MCCollection::MC12c) 
			JES_config_file = "../../ApplyJetCalibration/data/CalibrationConfigs/JES_Full2012dataset_Preliminary_MC12b_Sep23.config";

    	calibrationJES = new JetAnalysisCalib::JetCalibrationTool(jetAlgo,JES_config_file, !isMC);
  	}

}
void Correction::CalibrateJet(D3PDReader::JetD3PDObject * jet, Int_t dataYear, Double_t rhoKt4EM, Double_t mu, Double_t NPV)
{
	// just to Print a message
	if(printInfoJetCal == 0)
	{
		cout<<"--------------------------------------"<<endl;
		cout<<"Jet Calibration Applied"<<endl;
		cout<<"--------------------------------------"<<endl;		
		printInfoJetCal++;
	}

	for(Int_t i = 0; i < jet->n(); i++)
	{
		// Getting jet info
		Double_t Eraw    = 0.0;
		Double_t Ptraw   = 0.0;
		Double_t eta     = 0.0;
		Double_t eta_det = 0.0;		
		Double_t phi     = 0.0;
		Double_t m       = 0.0;
		Double_t Ax      = 0.0;
		Double_t Ay      = 0.0;
		Double_t Az      = 0.0;
		Double_t Ae      = 0.0;
		Double_t rho     = 0.0;
		if (dataYear==2012){
			Eraw    = (*jet)[i].emscale_E();
			Ptraw   = (*jet)[i].emscale_pt();
			eta     = (*jet)[i].emscale_eta();
			phi     = (*jet)[i].emscale_phi();
			m       = (*jet)[i].emscale_m();
	
			Ax 		= (*jet)[i].ActiveAreaPx();
			Ay      = (*jet)[i].ActiveAreaPy();
			Az      = (*jet)[i].ActiveAreaPz();
			Ae      = (*jet)[i].ActiveAreaE();
			rho     = rhoKt4EM;
		}
		else if(dataYear == 2011)
		{
			Eraw    = (*jet)[i].emscale_E();
			eta_det = (*jet)[i].emscale_eta();
			eta     = (*jet)[i].EtaOrigin();
			phi     = (*jet)[i].PhiOrigin();
			m       = (*jet)[i].MOrigin();
		}
		
		Double_t Ecal   = Eraw;
		Double_t pTcal  = Ptraw;
		Double_t etacal = eta;
		Double_t phical = phi;
		if(jetCalibration)
		{
			TLorentzVector jetCal;
			if(dataYear == 2012) jetCal = calibrationJES->ApplyJetAreaOffsetEtaJES(Eraw,eta,phi,m,Ax,Ay,Az,Ae,rho,mu,NPV);
			else if (dataYear == 2011) jetCal = calibrationJES->ApplyOffsetEtaJES(Eraw,eta_det,eta,phi,m,mu,NPV);

			pTcal  = jetCal.Pt();
    		etacal = jetCal.Eta();
    		phical = jetCal.Phi();
    		Ecal = jetCal.E();

			//if(isDebugCall)
			//{
			//	cout<<"Jets before smearing: "<< i;
			//	cout<<" Eta: "<<eta;
			//	cout<<" Phi: "<<phi;
			//	cout<<" Pt: "<<(*jet)[i].pt();
			//	cout<<" m: "<<m;
			//	cout<<" mu: "<<mu;
			//	cout<<" NPV: "<<NPV<<endl;

			//	cout<<"Other Inputs calibration function";
			//	cout<<" Eraw "<< Eraw;
			//	cout<<" eta_det "<< eta_det<<endl;;
			//	cout<<" Eraw "<< Eraw;
			//	cout<<" Eraw "<< Eraw;
			//	cout<<" Ax "<< Ax;
			//	cout<<" Ay "<< Ay;
			//	cout<<" Az "<< Az;
			//	cout<<" Ae "<< Ae;
			//	cout<<" rho "<< rho<<endl;
			//}


			(*jet)[i].E() = Ecal;
    		(*jet)[i].pt() = pTcal;
    		(*jet)[i].eta() = etacal;
    		(*jet)[i].phi() = phical;
			// To compare with Valerio
    		(*jet)[i].m() = jetCal.M();

			//if(isDebugCall)
			//{
			//	cout<<"Jets after smearing: "<< i;
			//	cout<<" Eta: "<<(*jet)[i].eta();
			//	cout<<" Phi: "<<(*jet)[i].phi();
			//	cout<<" Pt: "<<(*jet)[i].pt();
			//	cout<<" m: "<<(*jet)[i].m();
			//	cout<<" E: "<<(*jet)[i].E();
			//	cout<<endl;
			//}

			
		}
	}
}

////////////////////////////////////////////////////////////////////////////////////////
//							D0/Z0 smearing
////////////////////////////////////////////////////////////////////////////////////////
//	Init
void Correction::InitD0Z0Smear()
{
	TFile *d0_smearing_file = new TFile("InputFile/SmearingD0/impact_parameter_smearing.root");
    std::cout << "Setting up d0 smearing config file"<<std::endl;
    smearD0[0] = (TH2F*)d0_smearing_file->Get("smearD0_0");
    smearD0[1] = (TH2F*)d0_smearing_file->Get("smearD0_1");
    smearD0[2] = (TH2F*)d0_smearing_file->Get("smearD0_2");
    smearD0[0]->SetDirectory(0);
    smearD0[1]->SetDirectory(0);
    smearD0[2]->SetDirectory(0);
    smearD0X = smearD0[0]->GetXaxis();
    smearD0Y = smearD0[0]->GetYaxis();
    d0_smearing_file->Close();

    TFile *z0_smearing_file = new TFile("InputFile/SmearingD0/impact_parameter_smearing.root");
    std::cout << "Setting up z0 smearing confi file"<<std::endl;
    smearZ0[0] = (TH2F*)z0_smearing_file->Get("smearZ0_0");
    smearZ0[1] = (TH2F*)z0_smearing_file->Get("smearZ0_1");
    smearZ0[2] = (TH2F*)z0_smearing_file->Get("smearZ0_2");
    smearZ0[0]->SetDirectory(0);
    smearZ0[1]->SetDirectory(0);
    smearZ0[2]->SetDirectory(0);
    smearZ0X = smearZ0[0]->GetXaxis();
    smearZ0Y = smearZ0[0]->GetYaxis();
    z0_smearing_file->Close();  
}

// D0/Z0 Smearing
// For Electron
void Correction::SmearD0Z0(D3PDReader::ElectronD3PDObject * el,Int_t EventNumber)
{
	// just to Print a message
	if(printInfoD0Z0Smear == 0)
	{
		cout<<"--------------------------------------"<<endl;		
		cout<<"D0Z0 Smearing Applied (El type)"<<endl;
		cout<<"--------------------------------------"<<endl;		
		printInfoD0Z0Smear ++;
	}
	// Electron smearing
	for(Int_t i = 0; i < el->n(); i++)
	{
		Double_t el_d0   =(*el)[i].trackd0pvunbiased();
    	Double_t el_z0   =(*el)[i].trackz0pvunbiased();
    	Int_t nBL        =(*el)[i].nBLHits();
    	Double_t eta     =(*el)[i].tracketa();
    	Double_t pt      =(*el)[i].trackpt();

		// Fixing the values
		(*el)[i].trackd0pvunbiased() = el_d0 - 2.e-3 + GetD0SmearSigma(EventNumber, nBL, pt, eta, i);
		(*el)[i].trackz0pvunbiased() = el_z0+ GetZ0SmearSigma(EventNumber, nBL, pt, eta, i);
		//cout<<"After Smear: "<<(*el)[i].trackz0pvunbiased()<<" Before Smear: "<< el_z0 <<endl;
	}
}
// For Muon
void Correction::SmearD0Z0(D3PDReader::MuonD3PDObject * mu, Int_t EventNumber, Int_t type)
{
	// just to Print a message
	if(printInfoD0Z0Smear == 0)
	{
		cout<<"--------------------------------------"<<endl;		
		cout<<"D0Z0 Smearing Applied (Mu type)"<<endl;
		cout<<"--------------------------------------"<<endl;		
		printInfoD0Z0Smear ++;
	}

	for(Int_t i = 0; i < mu->n(); i++)
	{
		Double_t mu_d0   =(*mu)[i].trackd0pvunbiased();
    	Double_t mu_z0   =(*mu)[i].trackz0pvunbiased();
    	Int_t nBL        =(*mu)[i].nBLHits();
    	Double_t eta     = 0;
    	Double_t pt      = 0;
		if(type == leptonType::MuonCalo)
		{
			Double_t theta = (*mu)[i].id_theta();
			eta = -log(tan(theta/2));
			pt = (1/fabs((*mu)[i].id_qoverp())*sin((*mu)[i].id_theta()));
		}
		else if(type == leptonType::MuonStaco)
		{
			if((*mu)[i].isStandAloneMuon())
			{
				eta = (*mu)[i].eta();
				pt = (*mu)[i].pt();
			}
			else
			{
				Double_t theta = (*mu)[i].id_theta();
				eta = -log(tan(theta/2));
				pt = (1/fabs((*mu)[i].id_qoverp())*sin((*mu)[i].id_theta()));
			}
		}
		else {cout<<"Muon D0Z0 smear: Type not recognized";}
		
		// Fixing the values Do nothing for standalone muon
		if(!((*mu)[i].isStandAloneMuon()))
		{
			(*mu)[i].trackd0pvunbiased() = mu_d0 - 2.e-3 + GetD0SmearSigma(EventNumber, nBL, pt, eta, i);
			(*mu)[i].trackz0pvunbiased() = mu_z0 + GetZ0SmearSigma(EventNumber, nBL, pt, eta, i);
		}

		//if(isDebugCall)
		//{
		//	cout<<"----------------------"<<endl;
		//	cout<<"Muon "<<i<<endl;
		//	cout<<"Before mu_d0: "<<mu_d0<<endl;
		//	cout<<"Before mu_z0: "<<mu_z0<<endl;
		//	cout<<"id_qoverp: "<<(*mu)[i].id_qoverp()<<endl;
		//	cout<<"id_theta: "<<(*mu)[i].id_theta()<<endl;
		//	if(type == leptonType::MuonCalo) cout<<"Type: Calo Muon"<<endl;
		//	else if(type == leptonType::MuonStaco && !(*mu)[i].isStandAloneMuon()) cout<<"Type: Staco Muon"<<endl;
		//	else if(type == leptonType::MuonStaco && (*mu)[i].isStandAloneMuon() ) cout<<"Type: StandAlone Muon"<<endl;
		//	cout<<"nBL: "<< nBL<<endl;
		//	cout<<"eta: "<<eta<<endl;
		//	cout<<"pt: "<<pt<<endl;
		//	cout<<"After mu_d0: "<<(*mu)[i].trackd0pvunbiased()<<endl;
		//	cout<<"After mu_z0: "<<(*mu)[i].trackz0pvunbiased()<<endl;
		//}

	}
}

//D0 SmearSigma 
Double_t Correction::GetD0SmearSigma(Int_t EventNumber, Int_t nBL, Double_t pt, Double_t eta, Int_t index)
{
	smearD0Rand.SetSeed(EventNumber + 100*index);
	Double_t smearD0Val = 0;
	if(nBL >= 0)
	{
		if(nBL >= 2) nBL = 2;
		Double_t sinTheta = 1./cosh(eta);
	//	Double_t p = pt*cosh(eta);
		Double_t p_quant = 1./sqrt(pt*pt*sinTheta)/1000.;
		Int_t xBin = smearD0X->FindFixBin(eta);
		Int_t yBin = smearD0Y->FindFixBin(p_quant);
		Double_t sigma = smearD0[nBL]->GetBinContent(xBin, yBin);

		smearD0Val = smearD0Rand.Gaus(0, sigma);
	}
	return smearD0Val;
}

//Z0 SmearSigma 
Double_t Correction::GetZ0SmearSigma(Int_t EventNumber, Int_t nBL, Double_t pt, Double_t eta, Int_t index)
{
	smearZ0Rand.SetSeed(EventNumber + 100*index);
	Double_t smearZ0Val = 0;
	if(nBL >= 0)
	{
		if(nBL >= 2) nBL = 2;
		Double_t sinTheta = 1./cosh(eta);
	//	Double_t p = pt*cosh(eta);
		Double_t p_quant = 1./sqrt(pt*pt*sinTheta)/1000.;
		Int_t xBin = smearZ0X->FindFixBin(eta);
		Int_t yBin = smearZ0Y->FindFixBin(p_quant);
		Double_t sigma = smearZ0[nBL]->GetBinContent(xBin, yBin);

		smearZ0Val = smearZ0Rand.Gaus(0, sigma);
	}
	return smearZ0Val;
}

void Correction::ClearVars()
{
	bfEP_cl_Et.clear();
	bfEP_cl_pt.clear();
	electronEpErr.clear();
}
