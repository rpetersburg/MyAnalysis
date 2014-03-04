#include <stdlib.h>
#include <string>
#include "MyAnalysis/MuonObject.h"
#include <iostream>


using namespace std;
////////////////////////////////////////////////////////////////////////////////////////
//				Constructor and Destructor
////////////////////////////////////////////////////////////////////////////////////////
// Just The constructor
MuonObject::MuonObject(Int_t year, Bool_t tLooseCutZ4l)
{
	dataYear = year;
	looseCutZ4l = tLooseCutZ4l;

	// Init for consitency
	Hist = 0;

	// Output to help people
	cbAuthor = 0;
	saAuthor = 0;
	caloAuthor = 0;
}

MuonObject::~MuonObject()
{

}

// Fill the initial vector with muon objects, Sets the type as well
void MuonObject::FillMuon(D3PDReader::MuonD3PDObject * mu_branch, Int_t type, vector<Double_t> muEff, Bool_t isMC)
{
	for(Int_t i = 0; i < mu_branch->n(); i++)
	{
		ChargedLepton *temp = new ChargedLepton (&((*mu_branch)[i]), type, muEff[i], i, isMC);
		muInitEvent.push_back(temp);
	}
	// sanity check
	if(mu_branch->n() != muEff.size())
	{
		cout<<"--------------------------------------"<<endl;
		cout<<"ERROR: MuonObject::FillMuon; Eff vector and muon vector not the same size"<<endl;
		cout<<"Type: "<<type<<" Muon Vector: "<<muInitEvent.size()<<" Mu eff: "<<muEff.size()<<endl;		
		cout<<"--------------------------------------"<<endl;
	}
}
// Performs the Muon cut by calling the right function for the type
Bool_t MuonObject::MuonCut(Int_t *cutMuPass)
{
	muBfOverlap.clear();

	// Check for the Hist
	if(Hist == 0) 
	{
		cout<<"Muon Histogram not Set";
		return false;
	}
	for(vector<ChargedLepton *>::iterator itr = muInitEvent.begin();
			itr != muInitEvent.end(); ++itr)
	{
		ChargedLepton *mu_lep_i = *itr;
		D3PDReader::MuonD3PDObjectElement *mu_i = mu_lep_i->GetMuon();
		
		// Variable for cosmic muon cuts
		Double_t d0Cut = 1.;
		Double_t z0Cut = 10.;
		//Staco
		if(mu_lep_i->type == leptonType::MuonStaco)
		{
			if(!CutStaco(mu_i, cutMuPass, d0Cut, z0Cut)) continue;
		}
		//StandAlone
		else if(mu_lep_i->type == leptonType::MuonStandAlone)
		{
			if(!CutStandAlone(mu_i, cutMuPass, d0Cut, z0Cut)) continue;
			
		}
		//Calo
		else if(mu_lep_i->type == leptonType::MuonCalo)
		{
			if(!CutCalo(mu_i, cutMuPass, d0Cut, z0Cut)) continue;
		}
		// For Storing the events that have passed the cuts but not the Overlap
		// Temp Storage
		muBfOverlap.push_back(mu_lep_i);
	}

	// To Remove Overlap 
	RemoveOverlap(cutMuPass);
	if(muOverlapGoodEvent.size() > 0) return true;
	else return false;
}

// Currently only removes mu-mu overlap
// Other overlaps are implemented outside this class
// May 31th, 2013 3:00pm
void MuonObject::RemoveOverlap(Int_t * cutMuPass)
{
	muOverlapGoodEvent.clear();

	for(vector<ChargedLepton *>::iterator itr = muBfOverlap.begin();
			itr != muBfOverlap.end(); ++itr)
	{
		ChargedLepton *mu_lep_i = *itr;
		D3PDReader::MuonD3PDObjectElement *mu_i = mu_lep_i->GetMuon();
		Bool_t reject = false;
		
		// Calo
		if(mu_lep_i->type == leptonType::MuonCalo)
		{
			reject = OverlapCalo(mu_i);
			if(reject) caloAuthor++;
		}
		// Stand Alone
		else if(mu_lep_i->type == leptonType::MuonStandAlone)
		{
			reject = OverlapStandAlone(mu_i);
			if(reject) saAuthor++;
			
		}
		
		// If any are rejects, doesn't push the muon on the next Vector
		if(!reject)
		{
		//	cutMuPass[cutMuFlow::OverLap]++;
			muOverlapGoodEvent.push_back(mu_lep_i);
		}
	}
}

// To Check Overlap with Calo
Bool_t MuonObject::OverlapCalo(D3PDReader::MuonD3PDObjectElement *mu_curr)
{
	for(vector<ChargedLepton *>::iterator itr_j = muBfOverlap.begin();
			itr_j != muBfOverlap.end(); ++itr_j)
	{
		ChargedLepton *mu_lep_j = *itr_j;
		D3PDReader::MuonD3PDObjectElement *mu_j = mu_lep_j->GetMuon();
		if(mu_curr == mu_j) continue;
	
		// compare it to any Staco (ST+CB and StandAlone)
		if(mu_lep_j->type != leptonType::MuonCalo )
		{
			Double_t currEta = -TMath::Log(TMath::Tan(mu_curr->id_theta()*0.5));
			Double_t currPhi = mu_curr->id_phi();

			Double_t consEta = -TMath::Log(TMath::Tan(mu_j->id_theta()*0.5));
			Double_t consPhi = mu_j->id_phi();
			Hist->muOverlapHist[leptonType::MuonCalo]->Fill(DeltaR(currEta, currPhi, consEta, consPhi), Hist->weight);
			if(DeltaR(currEta, currPhi, consEta, consPhi) < 0.1) 
				return true;	
		}

	}

	return false;
}

// To Check Overlap with StandAlone
Bool_t MuonObject::OverlapStandAlone(D3PDReader::MuonD3PDObjectElement *mu_curr)
{
	for(vector<ChargedLepton *>::iterator itr_j= muBfOverlap.begin();
			itr_j != muBfOverlap.end(); ++itr_j)
	{
		ChargedLepton *mu_lep_j = *itr_j;
		D3PDReader::MuonD3PDObjectElement *mu_j = mu_lep_j->GetMuon();
		if(mu_curr == mu_j) continue;

		// Compare it to segment tagged muon
		if(mu_lep_j->type != leptonType::MuonCalo && mu_j->isSegmentTaggedMuon() == 1)
		{
			Double_t currEta = mu_curr->eta();
			Double_t currPhi = mu_curr->phi();

			Double_t consEta = mu_j->eta();
			Double_t consPhi = mu_j->phi();

			Double_t deltaRData = DeltaR(currEta, currPhi, consEta, consPhi);
			Hist->muOverlapHist[leptonType::MuonStandAlone]->Fill(deltaRData, Hist->weight);
				
			if(deltaRData < 0.2) 
				return true;	
		}

	}

	return false;
}

/////////////////////////////////////////////////////////////////////////////////
//						Helper Funtions
/////////////////////////////////////////////////////////////////////////////////
// Just to calculate DeltaR
Double_t MuonObject::DeltaR (Double_t eta_1, Double_t phi_1, Double_t eta_2, Double_t phi_2)
{
	Double_t dR=0;
	Double_t eta2 = (eta_1-eta_2)*(eta_1-eta_2);
	Double_t tmp_dphi = (fabs(phi_1-phi_2) > TMath::Pi()) ? 2*TMath::Pi()-fabs(phi_1-phi_2) : fabs(phi_1-phi_2);
	Double_t phi2 = tmp_dphi*tmp_dphi;
	dR = sqrt( eta2 + phi2 );
	return dR;
}

// Staco Muons cuts
Bool_t MuonObject::CutStaco(D3PDReader::MuonD3PDObjectElement *mu_i, Int_t *cutMuPass, 
		Double_t d0Cut, Double_t z0Cut)
{
	//Author cut
	Hist->muAuthorHist[leptonType::MuonStaco]->Fill(mu_i->author(), Hist->weight);
	if((mu_i->author() == 6 || mu_i->author() == 7) && mu_i->isStandAloneMuon() == 0)
	{cutMuPass[cutMuFlow::Author]++;}
	else return false;

	// Kinematic Cut 
	// Pt
	Double_t pTCut = 6000;
	if(looseCutZ4l) pTCut = 4000;
	Hist->muPTHist[leptonType::MuonStaco]->Fill(mu_i->pt(), Hist->weight); 
	if(mu_i->pt() > pTCut) 
	{cutMuPass[cutMuFlow::Pt]++;}
	else return false;

	// Eta
	Hist->muEtaHist[leptonType::MuonStaco]->Fill(mu_i->eta(), Hist->weight);	 
	if(fabs(mu_i->eta()) <  2.7)
	{cutMuPass[cutMuFlow::Eta]++;}
	else return false;
	
	//ID Cuts
	Int_t nBLayer = 0;
	Int_t nPix = 0;
	Int_t nSCT = 0;
	Int_t nHoles = 0;
	Double_t nLwrEta = 0;
	Double_t nUprEta = 0;
	if(dataYear == 2011)
	{
		nBLayer = 0;
		nPix = 1;
		nSCT = 5;
		nHoles = 3;
		nLwrEta = 0;
		nUprEta = 1.9;
	}
	else if (dataYear == 2012)
	{
		nBLayer = 0;
		nPix = 0;
		nSCT = 4;
		nHoles = 3;
		nLwrEta = 0.1;
		nUprEta = 1.9;
	}
	// Blayer
	if(dataYear == 2012 || (dataYear == 2011 && (!mu_i->expectBLayerHit() || mu_i->nBLHits() > nBLayer)))
	{cutMuPass[cutMuFlow::BLayer]++;}
	else return false;
	//Pix
	if((mu_i->nPixHits() + mu_i->nPixelDeadSensors ()) > nPix)
	{cutMuPass[cutMuFlow::Pix]++;}
	else return false;
	//STC
	if((mu_i->nSCTHits() + mu_i->nSCTDeadSensors ()) > nSCT)
	{cutMuPass[cutMuFlow::SCT]++;}
	else return false;
	// Holes
	if((mu_i->nPixHoles() + mu_i->nSCTHoles()) < nHoles)
	{cutMuPass[cutMuFlow::Holes]++;}
	else return false;
	// TRT
	Int_t nTRT = mu_i->nTRTHits() + mu_i->nTRTOutliers();
	Double_t dataEta = fabs(mu_i->eta());
	if((dataEta > nLwrEta && dataEta < nUprEta))
	{
		if(nTRT > 5 && mu_i->nTRTOutliers() < nTRT*0.9){}
		else return false;
	}
	else if(dataYear == 2011 && (dataEta < nLwrEta || dataEta >= nUprEta))
	{
		if(nTRT > 5)
		{
			if(mu_i->nTRTOutliers() < nTRT*0.9){}
			else return false;
		}
	}
	cutMuPass[cutMuFlow::TRT]++;
	
	//d0 and z0 cut
	Hist->muD0Hist[leptonType::MuonStaco]->Fill(mu_i->trackd0pvunbiased(), Hist->weight);
	Hist->muZ0Hist[leptonType::MuonStaco]->Fill(mu_i->trackz0pvunbiased(), Hist->weight);
	if(fabs(mu_i->trackd0pvunbiased()) < d0Cut &&
			fabs(mu_i->trackz0pvunbiased()) < z0Cut)
	{cutMuPass[cutMuFlow::D0]++;}
	else return false;
	
	return true;
}

// Calo muon cuts
Bool_t MuonObject::CutCalo(D3PDReader::MuonD3PDObjectElement *mu_i, Int_t *cutMuPass,
		Double_t d0Cut, Double_t z0Cut)
{
	// Author Cut
	Hist->muAuthorHist[leptonType::MuonCalo]->Fill(mu_i->author(), Hist->weight);	
	if(mu_i->author() == 16 && (mu_i->caloMuonIdTag() > 10 || mu_i->caloLRLikelihood() > 0.9))
	{cutMuPass[cutMuFlow::Author]++;}
	else return false;

	// Kinematic Cut
	Hist->muPTHist[leptonType::MuonCalo]->Fill(mu_i->pt(), Hist->weight);
	if(mu_i->pt() > 15000)
	{cutMuPass[cutMuFlow::Pt]++;}
	else return false;
	
	// Eta
	Hist->muEtaHist[leptonType::MuonCalo]->Fill(mu_i->eta(), Hist->weight);	
	if(fabs(mu_i->eta()) < 0.1)
	{cutMuPass[cutMuFlow::Eta]++;}
	else return false;

	//ID Cuts
	Int_t nBLayer = 0;
	Int_t nPix = 0;
	Int_t nSCT = 0;
	Int_t nHoles = 0;
	Double_t nEta = 0;
	if(dataYear == 2011)
	{
		nBLayer = 0;
		nPix = 1;
		nSCT = 5;
		nHoles = 3;
		nEta = 0.1;
	}
	else if (dataYear == 2012)
	{
		nBLayer = 0;
		nPix = 0;
		nSCT = 4;
		nHoles = 3;
		nEta = 0.1;
	}
	// Blayer
	if(dataYear == 2012 || (dataYear == 2011 && (!mu_i->expectBLayerHit() || mu_i->nBLHits() > nBLayer)))
	{cutMuPass[cutMuFlow::BLayer]++;}
	else return false;
	//Pix
	if((mu_i->nPixHits() + mu_i->nPixelDeadSensors ()) > nPix)
	{cutMuPass[cutMuFlow::Pix]++;}
	else return false;
	//STC
	if((mu_i->nSCTHits() + mu_i->nSCTDeadSensors ()) > nSCT)
	{cutMuPass[cutMuFlow::SCT]++;}
	else return false;
	// Holes
	if((mu_i->nPixHoles() + mu_i->nSCTHoles()) < nHoles)
	{cutMuPass[cutMuFlow::Holes]++;}
	else return false;
	//TRT
	Int_t nTRT = mu_i->nTRTHits() + mu_i->nTRTOutliers();
	Double_t dataEta = fabs(mu_i->eta());
	if(dataYear == 2011 && dataEta < 0.1)
	{
		if(nTRT < 6 || mu_i->nTRTOutliers() < nTRT *0.9){}
		else return false;
	}
	cutMuPass[cutMuFlow::TRT]++;

	//d0 and z0 cut
	Hist->muD0Hist[leptonType::MuonCalo]->Fill(mu_i->trackd0pvunbiased(), Hist->weight);
	Hist->muZ0Hist[leptonType::MuonCalo]->Fill(mu_i->trackz0pvunbiased(), Hist->weight);
	if(fabs(mu_i->trackd0pvunbiased()) < d0Cut &&
			fabs(mu_i->trackz0pvunbiased()) < z0Cut)
	{cutMuPass[cutMuFlow::D0]++;}
	else return false;

	return true;
}
// StandAlone cuts
Bool_t MuonObject::CutStandAlone(D3PDReader::MuonD3PDObjectElement *mu_i, Int_t *cutMuPass,
		Double_t d0Cut, Double_t z0Cut)
{
	Hist->muAuthorHist[leptonType::MuonStandAlone]->Fill(mu_i->author(), Hist->weight);		
	if(mu_i->author() == 6  && mu_i->isStandAloneMuon() == 1)
	{cutMuPass[cutMuFlow::Author]++;}
	else return false;

	// Kinematic Cut
	Double_t pTCut = 6000;
	if(looseCutZ4l) pTCut = 4000;
	Hist->muPTHist[leptonType::MuonStandAlone]->Fill(mu_i->pt(), Hist->weight);
	if(mu_i->pt() > pTCut)
	{cutMuPass[cutMuFlow::Pt]++;}
	else return false;

	// Eta
	Hist->muEtaHist[leptonType::MuonStandAlone]->Fill(mu_i->eta(), Hist->weight);	
	if(fabs(mu_i->eta()) <  2.7 && fabs(mu_i->eta()) > 2.5)
	{cutMuPass[cutMuFlow::Eta]++;}
	else return false;
	
	//ID Cut
	if((mu_i->nCSCEtaHits()+mu_i->nCSCPhiHits()) > 0 
			&& mu_i->nMDTEMHits() > 0 && mu_i->nMDTEOHits() > 0){}
	else return false;

	// No corresponding cuts for StandAlone
	cutMuPass[cutMuFlow::BLayer]++;
	cutMuPass[cutMuFlow::Pix]++;
	cutMuPass[cutMuFlow::SCT]++;
	cutMuPass[cutMuFlow::Holes]++;
	cutMuPass[cutMuFlow::TRT]++;
	cutMuPass[cutMuFlow::D0]++;

	return true;

}
// Set the histrograms
void  MuonObject::SetHist (HistContainer *curr_Hist)
{
	Hist = curr_Hist;	
}
// Clear vars
void MuonObject::clearVars()
{
	while(!muInitEvent.empty()) delete muInitEvent.back(), muInitEvent.pop_back();

	muInitEvent.clear(); 
	muBfOverlap.clear();	
	muOverlapGoodEvent.clear(); 
	muEvent.clear();

}
