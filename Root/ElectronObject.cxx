#include <stdlib.h>
#include <string>
#include "MyAnalysis/ElectronObject.h"

#include <iostream>


using namespace std;
////////////////////////////////////////////////////////////////////////////////////////
//				Constructor and Destructor
////////////////////////////////////////////////////////////////////////////////////////
// Just The constructor
ElectronObject::ElectronObject(Int_t year, Bool_t tuseLikelihood)
{
	dataYear = year;
	useLikelihood = tuseLikelihood;

	if(dataYear == 2011)
	{
		//elID2011= new H4l2011Defs();

		//	TPython::Exec("import os");
		//	TPython::LoadMacro("../../ElectronPhotonSelectorTools/python/HelloWorld.py");
		//	TPython::Eval("ConfiguredTElectronIsEMSelector(ROOT.egammaPID.ElectronIDLoosePP,electronPIDmenu.menuH4l2011)"); 

		TString configH4l = "../../ElectronPhotonSelectorTools/python/ConfiguredTElectronIsEMSelectors.py";
		TPython::LoadMacro(configH4l.Data());

		elID2011new = static_cast<Root::TElectronIsEMSelector*>((void*)TPython::Eval("ConfiguredTElectronIsEMSelector(ROOT.egammaPID.ElectronIDLoosePP,electronPIDmenu.menuH4l2011)")); 
		elID2011new->initialize();
	}
	else if(dataYear == 2012  && !useLikelihood)
	{
		ml_2013 =new Root::TElectronMultiLeptonSelector(); 
		ml_2013->initialize();
	}
	else if(dataYear==2012 && useLikelihood)
	{
		elID2012 = new Root::TElectronLikelihoodTool();
		elID2012->setPDFFileName("../../ElectronPhotonSelectorTools/data/ElectronLikelihoodPdfs.root"); // use this root file!
		elID2012->setOperatingPoint(LikeEnum::Loose);
		elID2012->initialize();

	}

	// Init for consistency
	Hist = 0;
}

// Fill the intial vector with elon objects, Sets the type as well
void ElectronObject::FillElectron(D3PDReader::ElectronD3PDObject * el_branch, Int_t type, vector<Double_t> telSmearVal, vector<Double_t> elEff, vector<Double_t> elRes, vector<Double_t> elCLPt, Bool_t isMC)
{
	for(Int_t i = 0; i < el_branch->n(); i++)
	{
		ChargedLepton *temp = new ChargedLepton (&((*el_branch)[i]), type, elEff[i], telSmearVal[i], i, isMC, elRes[i], elCLPt[i]);
		elInitEvent.push_back(temp);
	}
	elSmearVal = telSmearVal;

	// sanity check
	if(elSmearVal.size() != elInitEvent.size())
	{
		cout<<"--------------------------------------"<<endl;
		cout<<"ERROR: ElectronObject::FillElectron; Smear vector and electron vector not the same size"<<endl;
		cout<<"--------------------------------------"<<endl;
	}
	// sanity check
	if(elEff.size() != elInitEvent.size())
	{
		cout<<"--------------------------------------"<<endl;
		cout<<"ERROR: ElectronObject::FillElectron; Eff vector and electron vector not the same size"<<endl;
		cout<<"Electron Vector: "<<elInitEvent.size()<<" El eff: "<<elEff.size()<<endl;
		cout<<"--------------------------------------"<<endl;
	}
	// sanity check
	if(elRes.size() != elInitEvent.size())
	{
		cout<<"--------------------------------------"<<endl;
		cout<<"ERROR: ElectronObject::FillElectron; elRes vector and electron vector not the same size"<<endl;
		cout<<"Electron Vector: "<<elInitEvent.size()<<" El Res: "<<elRes.size()<<endl;
		cout<<"--------------------------------------"<<endl;
	}


}
// Perfroms the Electron cut by calling the right funtion for the type
Bool_t ElectronObject::ElectronCut(Int_t *cutElPass, Int_t Npv, Bool_t useRelaxedLoose)
{
	// Check if Hist are intilaized
	if(Hist == 0)
	{
		cout<<"Electron Histogram Not Init"<<endl;
		return false;
	}
	elBfOverlap.clear();
	// Simple counter
	Int_t i = -1;
	for(vector<ChargedLepton *>::iterator itr = elInitEvent.begin();
			itr != elInitEvent.end(); ++itr)
	{
		ChargedLepton *el_lep_i = *itr;
		D3PDReader::ElectronD3PDObjectElement *el_i = el_lep_i->GetElectron();

		i++;

		// Variable for cosmic elon cuts
		Double_t z0Cut = 10.;
		//Staco
		if(el_lep_i->type == leptonType::ElectronGSF)
		{
			//if(!CutGSF(el_i, cutElPass, z0Cut,elSmearVal[i], Npv)) continue;
			if(!CutGSF(el_i, cutElPass, z0Cut,el_lep_i->smearVal, Npv, useRelaxedLoose)) continue;
		}

		elBfOverlap.push_back(el_lep_i);
	}
	RemoveOverlap(cutElPass, Npv);
	if(elOverlapGoodEvent.size() > 0) return true;
	else return false;
}

// Currently only removes mu-mu overlap
// Other overlaps are implemented outside this class
// May 31th, 2013 3:00pm
void ElectronObject::RemoveOverlap(Int_t * cutElPass, Int_t NPV)
{
	elOverlapGoodEvent.clear();
	elBfCCOverlap.clear();
	for(vector<ChargedLepton *>::iterator itr = elBfOverlap.begin();
			itr != elBfOverlap.end(); ++itr)
	{
		ChargedLepton *el_lep_i = *itr;
		D3PDReader::ElectronD3PDObjectElement *el_i = el_lep_i->GetElectron();
		Bool_t reject = false;

		reject = OverlapEE(el_i, NPV);

		// Overwrite if not loose 
		if(! AuthorCut(el_i, el_lep_i->smearVal, NPV, false)) reject = false;
		// If any are rejects, doesn't push the electron on the next Vector
		if(!reject)
		{
			elBfCCOverlap.push_back(el_lep_i);
			cutElPass[cutElFlow::OverLapElEl]++;
		}
	}

	for(vector<ChargedLepton *>::iterator itr = elBfCCOverlap.begin();
			itr != elBfCCOverlap.end(); ++itr)
	{
		ChargedLepton *el_lep_i = *itr;
		D3PDReader::ElectronD3PDObjectElement *el_i = el_lep_i->GetElectron();
		Bool_t reject = false;

		reject = OverlapEECC(el_i, NPV);
		
		if(! AuthorCut(el_i, el_lep_i->smearVal, NPV, false)) reject = false;
		// If any are rejects, doesn't push the electron on the next Vector
		if(!reject)
		{
			elOverlapGoodEvent.push_back(el_lep_i);
			cutElPass[cutElFlow::OverLapClElEl]++;
		}
	}

}

// To Check Overlap with ee overlap 
Bool_t ElectronObject::OverlapEE(D3PDReader::ElectronD3PDObjectElement *el_curr, Int_t NPV)
{
	// Getting the Current Variables
	Double_t d0_curr = -1;
	Double_t z0_curr = -1;
	Double_t phi_curr = -1;
	Double_t qoverp_curr = -1;
	Double_t Eta_curr = -1;
	Double_t Et_curr = -1;

	if(dataYear == 2011)
	{
		d0_curr = el_curr->trackd0();
		z0_curr = el_curr->trackz0();
		phi_curr = el_curr->trackphi();
		qoverp_curr = el_curr->trackqoverp();
		Eta_curr = el_curr->tracketa();
	}
	else if (dataYear == 2012)
	{
		d0_curr = el_curr->Unrefittedtrack_d0();
		z0_curr = el_curr->Unrefittedtrack_z0();
		phi_curr = el_curr->Unrefittedtrack_phi();
		qoverp_curr = el_curr->Unrefittedtrack_qoverp();
		Eta_curr = el_curr->Unrefittedtrack_eta();
	}
	Et_curr = el_curr->cl_E()/cosh(Eta_curr);

	for(vector<ChargedLepton *>::iterator itr_j = elBfOverlap.begin();
			itr_j != elBfOverlap.end(); ++itr_j)
	{
		ChargedLepton *el_lep_j = *itr_j;
		D3PDReader::ElectronD3PDObjectElement *el_j = el_lep_j->GetElectron();
		if(el_curr == el_j) continue;

		// Continue if the current electron is not a loose electron
		if(! AuthorCut(el_j, el_lep_j->smearVal, NPV, false)) continue;
		
		/// For the electron to compare
		Double_t d0_j = -1;
		Double_t z0_j = -1;
		Double_t phi_j = -1;
		Double_t qoverp_j = -1;
		Double_t Eta_j = -1;
		Double_t Et_j = -1;

		if(dataYear == 2011)
		{
			d0_j = el_j->trackd0();
			z0_j = el_j->trackz0();
			phi_j = el_j->trackphi();
			qoverp_j = el_j->trackqoverp();
			Eta_j = el_j->tracketa();
		}
		else if (dataYear == 2012)
		{
			d0_j = el_j->Unrefittedtrack_d0();
			z0_j = el_j->Unrefittedtrack_z0();
			phi_j = el_j->Unrefittedtrack_phi();
			qoverp_j = el_j->Unrefittedtrack_qoverp();
			Eta_j = el_j->Unrefittedtrack_eta();
		}
		Et_j = el_j->cl_E()/cosh(Eta_j);

		// Comparing they share the same ID and if the curr has lower ET, reject it
		if( d0_curr == d0_j &&
				z0_curr == z0_j &&
				phi_curr == phi_j &&
				qoverp_curr == qoverp_j &&
				Et_curr < Et_j
		  ) return true;
	}

	return false;
}
// To check for ee cluster over lap
Bool_t ElectronObject::OverlapEECC(D3PDReader::ElectronD3PDObjectElement *el_curr, Int_t NPV)
{
	Double_t phi_curr = -1;
	Double_t eta_curr = -1;
	Double_t eta_track_curr = -1;
	Double_t Et_curr = -1;

	Double_t cutPhi = -1;
	Double_t cutEta = -1;

	if(dataYear == 2011)
	{
		cutPhi = 0;
		cutEta = 0;
	}
	else if(dataYear == 2012)
	{
		cutPhi = 5*0.025;
		cutEta = 3*0.025;
	}

	phi_curr = el_curr->cl_phi();
	eta_curr = el_curr->cl_eta();
	eta_track_curr = el_curr->tracketa();
	Et_curr = el_curr->cl_E()/cosh(eta_track_curr);

	for(vector<ChargedLepton *>::iterator itr_j = elBfOverlap.begin();
			itr_j != elBfOverlap.end(); ++itr_j)
	{
		ChargedLepton *el_lep_j = *itr_j;
		D3PDReader::ElectronD3PDObjectElement *el_j = el_lep_j->GetElectron();
		if(el_curr == el_j) continue;
		// Continue if the current electron is not a loose electron
		if(! AuthorCut(el_j, el_lep_j->smearVal, NPV, false)) continue;

		/// For the electron to compare
		Double_t phi_j = el_j->cl_phi();
		Double_t eta_j = el_j->cl_eta();
		Double_t eta_track_j = el_j->tracketa();
		Double_t Et_j = el_j->cl_E()/cosh(eta_track_j);

		Double_t deltaEta = fabs(eta_curr - eta_j);
		Double_t deltaPhi = fabs(phi_curr - phi_j);
		deltaPhi = (deltaPhi > TMath::Pi()) ? 2*TMath::Pi()-deltaPhi : deltaPhi;

		// Comparing they share the same ID and if the curr has lower ET, reject it
		if( deltaEta < cutEta &&
				deltaPhi < cutPhi &&
				Et_curr < Et_j
		  ) return true;
	}

	return false;

}


/////////////////////////////////////////////////////////////////////////////////
//						Helper Funtions
/////////////////////////////////////////////////////////////////////////////////
// Just to calculate DeltaR
Double_t ElectronObject::DeltaR (Double_t eta_1, Double_t phi_1, Double_t eta_2, Double_t phi_2)
{
	Double_t dR=0;
	Double_t eta2 = (eta_1-eta_2)*(eta_1-eta_2);
	Double_t tmp_dphi = (fabs(phi_1-phi_2) > TMath::Pi()) ? 2*TMath::Pi()-fabs(phi_1-phi_2) : fabs(phi_1-phi_2);
	Double_t phi2 = tmp_dphi*tmp_dphi;
	dR = sqrt( eta2 + phi2 );
	return dR;
}

// Staco Electrons cuts
Bool_t ElectronObject::CutGSF(D3PDReader::ElectronD3PDObjectElement *el_i, Int_t *cutElPass, Double_t z0Cut, Double_t smearVal_i, Int_t Npv, Bool_t useRelaxed)
{
	//Author cut
	Hist->elAuthorHist->Fill(el_i->author(), Hist->weight);
	if((el_i->author() == 1 || el_i->author() == 3))
	{cutElPass[cutElFlow::Author]++;}
	else return false;

	if(!AuthorCut(el_i, smearVal_i, Npv, useRelaxed)) return false;
	cutElPass[cutElFlow::Loose]++;

	//Kinematics
	Double_t Et=0.0;
	Double_t Eta_cl=0.0;
	Double_t Eta_tr=0.0;

	Eta_cl = el_i->cl_eta();
	Eta_tr = el_i->tracketa(); 
	Et = el_i->cl_E()/(cosh(Eta_tr));

	//Eta
	Hist->elEtaHist->Fill(Eta_cl, Hist->weight);
	if(fabs(Eta_cl) < 2.47) // These variables come from above intiia
	{cutElPass[cutElFlow::Eta]++;}
	else return false;
	// Et
	Hist->elETHist->Fill(Et, Hist->weight);
	if(Et > 7000)
	{cutElPass[cutElFlow::Et]++;}
	else return false;

	// Object Quality
	if((el_i->OQ() & 1446) == 0)
	{cutElPass[cutElFlow::ObjectQuality]++;}
	else return false;

	// Z0 Cut
	Hist->elZ0Hist->Fill(el_i->trackz0pvunbiased(), Hist->weight);
	if(fabs(el_i->trackz0pvunbiased()) < z0Cut)
	{cutElPass[cutElFlow::Z0]++;}
	else return false;	

	return true;
}

void ElectronObject::SetHist(HistContainer *curr_Hist)
{
	Hist = curr_Hist;	

}


Bool_t ElectronObject::AuthorCut(D3PDReader::ElectronD3PDObjectElement *el_i, Double_t smearVal_i, Int_t ip, Bool_t useRelaxed)
{
	Bool_t PassElectronID = true;

	Double_t Et=0.0;
	Double_t Et_cl=0.0;
	Double_t Et_cl_s2=0.0;
	Double_t pt=0.0;
	Double_t Eta_cl=0.0;
	Double_t Eta_cl_s2=0.0;
	Double_t Eta_trk=0.0;
	Int_t    Author=0;

	Double_t f3=0.0;
	Double_t rHad=0.0; 
	Double_t rHad1=0.0; 
	Double_t Had=0.0; 
	Double_t Had1=0.0;
	Double_t Reta=0.0; 
	Double_t w2=0.0;
	Double_t f1=0.0; 
	Double_t wstot=0.0; 
	Double_t DEmaxs1=0.0;
	Double_t deltaEta=0.0;
	Double_t deltaPhi=0.0;
	Double_t d0=0.0;
	Double_t TRratio=0.0; 
	Int_t    nSCT=0.0;
	Int_t    nSCTOutliers=0.0;
	Int_t    nTRThigh=0.0;
	Int_t    nTRThighOutliers=0.0;
	Int_t    nTRTXenonHits=0.0;
	Int_t    nTRT=0; 
	Int_t    nTRTOutliers=0;
	Int_t    nSi=0;
	Int_t    nSiOutliers=0;
	Int_t    nPix=0;
	Int_t    nPixOutliers=0; 
	Int_t    nBlayer=0;
	Int_t    nBlayerOutliers=0; 
	Bool_t   expectBlayer=false;
	Double_t rTRT=0.0;
	Double_t deltaPhiRescaled=0.0;
	Int_t    nSiDeadSensors=0;
	Int_t    nPixDeadSensors=0;
	Int_t    nBlayerHits=0;
	Double_t rphi=0.0;
	Int_t    convBit=0;
	Double_t d0sigma=0.0;
	Double_t ws3=0.0;
	Double_t e237=0.0;
	Double_t e277=0.0;
	Double_t w1=0.0;
	Double_t emax2=0.0;
	Double_t emax=0.0;
	Double_t emin=0.0;
	Double_t fracm=0.0;
	Double_t ep=0.0;

	Eta_cl = el_i->cl_eta();
	Eta_trk = el_i->tracketa(); 
	Eta_cl_s2 = el_i->etas2();
	Et = el_i->cl_E()/(cosh(Eta_trk));
	Et_cl = el_i->cl_E()/ (cosh(Eta_cl));
	Et_cl_s2= el_i->cl_E()/(cosh(Eta_cl_s2));

	Author = el_i->author();

	if(dataYear == 2011)
	{
		Et_cl_s2=Et_cl_s2/smearVal_i;
		//	rHad = el_i->Ethad()/Et_cl_s2; 
		//	rHad1 = el_i->Ethad1()/Et_cl_s2; 
		//	Reta = el_i->reta(); 
		//	w2 = el_i->weta2();
		//	f1 = el_i->f1(); 
		//	wstot=el_i->wstot(); 
		//	DEmaxs1=( el_i->emaxs1() - el_i->Emax2())/( el_i->emaxs1() + el_i->Emax2() );
		//	deltaEta=el_i->deltaeta1(); 
		//	nSi = el_i->nSiHits(); 
		//	nSiOutliers = el_i->nSCTOutliers()+el_i->nPixelOutliers();
		//	nPix = el_i->nPixHits () ; 
		//	nPixOutliers = el_i->nPixelOutliers();
		//	PassElectronID=elID2011->passH4l2011(Eta_cl_s2, 
		//    	                                 Et_cl_s2,
		//        	                             rHad, 
		//            	                         rHad1, 
		//                	                     Reta, 
		//                    	                 w2,
		//                        	             f1, 
		//                            	         wstot, 
		//                                	     DEmaxs1,
		//                                    	 deltaEta, 
		//                       	    	         nSi, 
		//                        	             nSiOutliers,
		//    	                                 nPix, 
		//        	                             nPixOutliers,
		//            	                         false,
		//                	                     false);
		e237=el_i->E237();
		e277=el_i->E277();
		Had=el_i->Ethad(); 
		Had1=el_i->Ethad1();
		w1=el_i->ws3();
		w2=el_i->weta2();
		f1=el_i->f1(); 
		emax2=el_i->Emax2();
		emax=el_i->emaxs1();
		emin=el_i->Emins1();
		wstot=el_i->wstot(); 
		fracm=el_i->fside(); 
		f3=el_i->f3(); 
		nBlayer=el_i->nBLHits();
		nBlayerOutliers=el_i->nBLayerOutliers();
		nPix=el_i->nPixHits(); 
		nPixOutliers=el_i->nPixelOutliers();
		nSCT=el_i->nSCTHits(); 
		nSCTOutliers=el_i->nSCTOutliers();
		nTRT=el_i->nTRTHits(); 
		nTRTOutliers=el_i->nTRTOutliers();
		nTRThigh= el_i->nTRTHighTHits();
		nTRThighOutliers= el_i->nTRTHighTOutliers();
		nTRTXenonHits=1.0;//el_i->nTRTXenonHits();    
		d0=el_i->trackd0_physics();
		deltaEta=el_i->deltaeta1(); 
		deltaPhi=el_i->deltaphi2(); 
		ep=1.;
		expectBlayer=el_i->expectHitInBLayer();
		EMAmbiguityType::AmbiguityResult amb=EMAmbiguityType::ELECTRON;

		rHad1=el_i->Ethad1()/Et_cl_s2;
		rHad=el_i->Ethad()/Et_cl_s2; 
		Reta=el_i->reta(); 
		//if (e277<0.) Reta=0.; // to be in sync with Athena macro
		//DEmaxs1=( el_i->emaxs1()-el_i->Emax2() )/ (el_i->emaxs1()+el_i->Emax2() );
		DEmaxs1= fabs(el_i->emaxs1()+el_i->Emax2())>0. ? ( el_i->emaxs1()-el_i->Emax2() )/( el_i->emaxs1()+el_i->Emax2() ): 0.;
		nSi=el_i->nSiHits(); 
		nSiOutliers=el_i->nSCTOutliers()+el_i->nPixelOutliers();

		Root::TAccept taccept;

		taccept=elID2011new->accept(Eta_cl_s2, 
				fabs(Et_cl_s2),
				e237,    
				e277,    
				Had1, 
				Had,
				w1,      
				w2,      
				f1, 
				emax2,    
				emax,    
				emin, 
				wstot,
				fracm,   
				f3, 
				nBlayer,
				nBlayerOutliers,
				nPix, 
				nPixOutliers,
				nSCT, 
				nSCTOutliers,
				nTRThigh,         
				nTRThighOutliers, 
				nTRT, 
				nTRTOutliers, 
				nTRTXenonHits,    
				d0,
				fabs(deltaEta), 
				deltaPhi, 
				ep,               
				expectBlayer,
				amb);
		PassElectronID=(bool)taccept;


	}
	else if(dataYear == 2012 && !useLikelihood)
	{
		// Uncorr  Et_cl_s2 for the isEM macro shower shape + Et bin
		Et_cl_s2=Et_cl_s2/smearVal_i;

		f3 = el_i->f3(); 
		rHad = el_i->Ethad()/Et_cl_s2; 
		rHad1 = el_i->Ethad1()/Et_cl_s2; 
		Reta = el_i->reta(); 
		w2 = el_i->weta2();
		f1 = el_i->f1(); 
		wstot = el_i->wstot(); 
		DEmaxs1 =  fabs(el_i->emaxs1()+el_i->Emax2())>0. ? ( el_i->emaxs1()-el_i->Emax2() )/( el_i->emaxs1() + el_i->Emax2() ): 0.;
		deltaEta = el_i->deltaeta1(); 
		d0 = el_i->trackd0_physics();
		TRratio = el_i->TRTHighTOutliersRatio(); 
		nTRT = el_i->nTRTHits(); 
		nTRTOutliers = el_i->nTRTOutliers();
		nSi = el_i->nSiHits(); 
		nSiOutliers = el_i->nSCTOutliers()+el_i->nPixelOutliers();
		nPix = el_i->nPixHits(); 
		// 
		nPixOutliers = el_i->nPixelOutliers();
		nBlayer = el_i->nBLHits();
		nBlayerOutliers = el_i->nBLayerOutliers(); 
		expectBlayer = el_i->expectHitInBLayer();  
		deltaPhiRescaled = el_i->deltaphiRescaled(); 
		nSiDeadSensors = el_i->nSCTDeadSensors() + el_i->nPixelDeadSensors(); 
		nPixDeadSensors = el_i->nPixelDeadSensors();
		nTRThigh = el_i->nTRTHighTHits();
		nTRThighOutliers = el_i->nTRTHighTOutliers();

		//nTRTOutliers = el_i->nTRTOutliers(); 
		nTRTOutliers = el_i->nTRTOutliers(); 
		nBlayerHits =  el_i->nBLHits();


		rTRT  =  (nTRT+nTRTOutliers) > 0 ?  ((double) (nTRThigh+nTRThighOutliers)/(nTRT+nTRTOutliers) ) : 0.;
		int nTRTTotal  =  nTRT+nTRTOutliers;

		double dpOverp  = 0;

		for (Int_t i =  0; i< (Int_t) el_i->refittedTrack_LMqoverp().size();++i){
			if( (el_i->refittedTrack_author().at(i) ==4)){
				dpOverp =  1- (el_i->trackqoverp()/( (el_i->refittedTrack_LMqoverp().at(i)) ));
			}
		}

		Root::TAccept taccept;

		taccept=ml_2013->accept(Eta_cl_s2, 
				Et_cl_s2,
				rHad, 
				rHad1, 
				Reta, 
				w2,
				f1,
				f3, 
				wstot, 
				DEmaxs1, 
				deltaEta, 
				nSi,
				nSiDeadSensors, 
				nPix, 
				nPixDeadSensors,
				deltaPhiRescaled,
				dpOverp,
				rTRT,
				nTRTTotal,
				nBlayerHits,
				expectBlayer,
				false);

		PassElectronID=(bool)taccept;
	}
	else if(dataYear == 2012 && useLikelihood)
	{
		Et_cl_s2=Et_cl_s2/smearVal_i;
		f3 = el_i->f3(); 
		rHad = el_i->Ethad()/Et_cl_s2; 
		rHad1 = el_i->Ethad1()/Et_cl_s2; 
		Reta = el_i->reta(); 
		w2 = el_i->weta2();
		w1 = el_i->weta2();
		f1 = el_i->f1(); 
		DEmaxs1 =  fabs(el_i->emaxs1()+el_i->Emax2())>0. ?  (el_i->emaxs1()-el_i->Emax2()) / (el_i->emaxs1()+el_i->Emax2()) : 0.;
		DEmaxs1 =  el_i->emaxs1()+el_i->Emax2() == 0. ? 0. :  (el_i->emaxs1()-el_i->Emax2()) / (el_i->emaxs1()+el_i->Emax2());

		deltaEta = el_i->deltaeta1(); 
		d0 = el_i->trackd0pvunbiased();
		// TRratio = el_i->TRTHighTHitsRatio(); 
		// Use nTRTHTOutliersRatio instead of HTHitsRatio.
		TRratio = el_i->TRTHighTOutliersRatio(); 
		d0sigma = el_i->tracksigd0pvunbiased();
		rphi = el_i->rphi();
		ws3 = el_i->ws3();
		deltaPhiRescaled =  el_i->deltaphiRescaled(); 
		nSi = el_i->nSiHits(); 
		//nSiOutliers = el_i->nSCTOutliers()+el_i->nPixelOutliers();
		// Use dead Si instead of outliers.
		nSiOutliers = el_i->nSCTDeadSensors()+el_i->nPixelDeadSensors();
		nPix = el_i->nPixHits(); 
		//nPixOutliers = el_i->nPixelOutliers();
		// Use dead pixels instead of outliers.
		nPixOutliers = el_i->nPixelDeadSensors();
		nBlayer = el_i->nBLHits();
		nBlayerOutliers = el_i->nBLayerOutliers();
		expectBlayer = el_i->expectHitInBLayer();
		if (el_i->expectHitInBLayer() == -999) expectBlayer = 1;
		//  egammaPID::ConversionMatch_Electron  =  = 1 
		convBit = el_i->isEM() & 0x1 << 1; 

		double dpOverp  = 0;
		if  (Author == 1 || Author == 3) {
			for (Int_t i =  0; i<el_i->refittedTrack_LMqoverp().size();++i){
				if (el_i->refittedTrack_author().at(i) == 4){
					dpOverp =  1- el_i->trackqoverp()/ el_i->refittedTrack_LMqoverp().at(i) ;
				}
			}
		}

	
		// For background control regions
		if(useRelaxed) elID2012->setOperatingPoint(LikeEnum::LooseRelaxed);
		else elID2012->setOperatingPoint(LikeEnum::Loose);

		Double_t discriminant = elID2012->calculate(Eta_cl_s2, 
				Et_cl_s2,
				f3, 
				rHad, 
				rHad1,
				Reta, 
				w2, 
				f1, 
				DEmaxs1,
				deltaEta, 
				d0, 
				TRratio,
				d0sigma,
				rphi,
				dpOverp ,
				deltaPhiRescaled,
				Double_t(ip));

		PassElectronID = elID2012->accept(discriminant,
				Eta_cl_s2, 
				Et_cl_s2,
				nSi,
				nSiOutliers, 
				nPix, 
				nPixOutliers,
				nBlayer, 
				nBlayerOutliers, 
				expectBlayer,
				convBit, 
				Double_t(ip));



	}
	return PassElectronID;
}


// To Clear the vars
void ElectronObject::clearVars()
{
	while(!elInitEvent.empty()) delete elInitEvent.back(), elInitEvent.pop_back();
	
	elInitEvent.clear(); 
	elBfOverlap.clear(); 		
	elBfCCOverlap.clear();
	elOverlapGoodEvent.clear(); 
	elEvent.clear();
	elSmearVal.clear();	
}
